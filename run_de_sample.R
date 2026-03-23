library(qs)
library(Seurat)
library(edgeR)
library(fdrtool)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(gridExtra)


seur <- qread("seur_multiome_final_harmony.qs")
seur <- subset(seur, assign_broad=='Singlet')
seur@meta.data = droplevels(seur@meta.data)
seur$CellType=factor(seur$CellType, levels=c('Car3','L2-3.IT','L4-5.IT','L5.ET','L5.IT','L5.NP','L6.CT','L6.IT','L6b',
	'Pvalb','Sst','Vip'))


# num singlets each guide across samples and cell types
tab = table(seur$sample, seur$assign, seur$CellType)
df_counts <- as.data.frame(tab)
colnames(df_counts) <- c("sample", "assignment", "cell_type", "count")
df_counts$cell_type <- as.factor(df_counts$cell_type)
df_counts$assignment <- as.factor(df_counts$assignment)

df_fractions <- df_counts %>% group_by(sample) %>% mutate(total_count = sum(count), fraction = count / total_count) %>% ungroup()

ggplot(df_fractions, aes(x = assignment, y = fraction)) +
  geom_bar(stat = "identity") +
  facet_grid(cell_type ~ sample) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    strip.text = element_text(size = 8), # Make facet labels readable
    legend.position = "none"             # Hide legend since X-axis labels indicate assignment
  ) +
  labs(
    x = "Assignment",
    y = "Fraction of cells"
  )
ggsave("barchart_frac_singlets_fix_y.pdf", height=12, width=16)

pb = AggregateExpression(seur,group.by = c("sample", "CellType", "assign"),
        return.seurat = FALSE,slot = "counts")$RNA
colnames(pb) = gsub("vivo_", "vivo", colnames(pb))
meta.pb <- data.frame(id = colnames(pb)) %>%
  tidyr::separate(id, into = c("sample", "CellType", "guide"), sep = "_", remove = FALSE)
rownames(meta.pb) = NULL
meta.pb <- meta.pb %>% tibble::column_to_rownames("id")
meta.pb$batch <- stringr::str_split(meta.pb$sample, "\\.", simplify=T)[,1]
meta.pb <- meta.pb %>% mutate(across(where(is.character), as.factor))
meta.pb$id = NULL
meta.pb$guide <- relevel(as.factor(meta.pb$guide), ref = "NT2-g2")

run_perturb_de <- function(cell_type_name, counts, metadata) {
  # subset
  idx <- metadata$CellType == cell_type_name
  sub_counts <- counts[, idx]
  sub_meta <- metadata[idx, ]
  
  # make sure we have more than one guide to compare
  if(length(unique(sub_meta$guide)) < 2) return(NULL)
  
  dge <- DGEList(counts = sub_counts, samples = sub_meta)
  keep <- filterByExpr(dge, group = sub_meta$guide)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  # model design
  design <- model.matrix(~ batch + guide, data = dge$samples)
  colnames(design) <- make.names(colnames(design))
  
  # fit
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  guide_cols <- grep("^guide", colnames(design), value = TRUE)
  
  # loop through each guide to compare against NT2-g2
  cell_type_results <- lapply(guide_cols, function(g_col) {
    res <- glmQLFTest(fit, coef = g_col)
    tab <- topTags(res, n = Inf)$table
    
    # fdrtool for empirical null correction
    ft <- fdrtool(tab$PValue, statistic = "pvalue", plot = FALSE)
    tab$fdrtool_q <- ft$qval
    
    # Metadata labels
    tab$gene <- rownames(tab)
    tab$cell_type <- cell_type_name
    tab$guide <- gsub("^guide", "", g_col)
    tab$comparison <- paste0(tab$guide, "_vs_NT2")
    
    return(tab)
  })
  
  return(do.call(rbind, cell_type_results))
}

all_cell_types <- levels(meta.pb$CellType)
res <- lapply(all_cell_types, run_perturb_de, counts = pb, metadata = meta.pb)
out <- do.call(rbind, res)

saveRDS(out, "table_DE_guide_vs_NT2-g2_v2.rds")
write.table(out, "table_DE_guide_vs_NT2-g2_v2.tsv", sep='\t', quote=F, row.names=F, col.names=T)

draw_celltype_volcano <- function(cell_type, comp_name, deg_data) {
        plot_df <- deg_data %>% filter(cell_type == !!cell_type, comparison == !!comp_name)
        if(nrow(plot_df) == 0) return(NULL)
        p <- EnhancedVolcano(plot_df,
                lab = plot_df$gene,
                x = 'logFC',
                y = 'fdrtool_q',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = paste(cell_type),
                subtitle = comp_name, # Display the drug dose here
                caption = paste(nrow(plot_df), "genes"),
                legendPosition = 'none',
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'))
        return(p)
}

plot_grid_info <- expand.grid(
        cell_type = unique(out$cell_type),
        comparison = unique(out$comparison),
        stringsAsFactors = F
) %>% arrange(cell_type)

for (ct in all_cell_types) {
cur_info = plot_grid_info %>% filter(cell_type == ct)
print(cur_info)
all_volcanos <- mapply(
        draw_celltype_volcano,
        cell_type = cur_info$cell_type,
        comp_name = cur_info$comparison,
        MoreArgs = list(deg_data = out),
        SIMPLIFY = FALSE
)

pdf(paste0("volcano_DE_guide_vs_NT2-g2_",ct,".pdf"), width=20, height=16)
print(do.call(grid.arrange, c(all_volcanos, ncol = 5)))
dev.off()
}

