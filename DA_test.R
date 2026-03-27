library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(tidyr)
library(glmGamPoi)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(patchwork)
library(qs)
library(glue)
library(matrixStats)
library(DESeq2)
source("~/kwanho/src/seurat_tools.R")

message("Start!")

seur <- qread("seur_multiome_final.qs")
seur <- subset(seur, assign %in% c('Doublet','Not assigned'), invert=T)

message("make pseudobulk!")

DefaultAssay(seur) <- "peaks"

seur$sample = gsub("_", "", seur$sample)
seur$bulk_id <- paste(seur$sample, seur$CellType, seur$assign, sep = "_")

pb_counts <- AggregateExpression(
  object = seur,
  group.by = "bulk_id",
  assays = "peaks",
  slot = "counts",
  return.seurat = FALSE
)$peaks

meta_pb <- seur@meta.data %>%
  dplyr::select(sample, CellType, assign) %>%
  distinct() %>%
  mutate(bulk_id = paste(sample, CellType, assign, sep = "_"))

rownames(meta_pb) = NULL
meta_pb = meta_pb %>% tibble::column_to_rownames("bulk_id")

rownames(meta_pb) <- gsub("_", "-", rownames(meta_pb))

pb_counts <- pb_counts[, rownames(meta_pb)]

##########################################################

run_da_deseq2 <- function(cell_type, counts, metadata) {
  require(DESeq2)
  # Subset
  sub_meta <- metadata[metadata$CellType == cell_type, ]
  sub_counts <- counts[, rownames(sub_meta)]

  # peak selection
  min_samples <- max(1, floor(ncol(sub_counts) * 0.05))
  keep_peaks <- rowSums(sub_counts > 0) >= min_samples
  sub_counts <- sub_counts[keep_peaks, ]
  
  sub_meta$assign <- relevel(as.factor(sub_meta$assign), ref = "NT2-g2")

  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(sub_counts),
    colData = sub_meta,
    design = ~ sample + assign
  )

  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)

  gRNAs <- levels(sub_meta$assign)
  targets <- gRNAs[gRNAs != "NT2-g2"]

  res_list <- lapply(targets, function(g) {
    res <- results(dds, contrast = c("assign", g, "NT2-g2"), tidy = TRUE)
    res$comparison <- g
    res$cell_type <- cell_type
    return(res)
  })
  
  return(do.call(rbind, res_list))
}

########################################################

message("DA test!")

cell_types <- unique(meta_pb$CellType)

da_results_list <- lapply(cell_types, function(ct) {
  
  # Safety Check: Count replicates per gRNA for this CellType
  sub_meta <- meta_pb[meta_pb$CellType == ct, ]
  
  if (nrow(sub_meta) <= length(unique(sub_meta$assign))) {
    message(paste("Skipping", ct, ": Not enough replicates for ~sample + assign design."))
    return(NULL)
  }

  message(paste("Running DA for:", ct))
  
  result <- tryCatch({
    run_da_deseq2(cell_type = ct, counts = pb_counts, metadata = meta_pb)
  }, error = function(e) {
    message(paste("Error in cell type", ct, ":", e$message))
    return(NULL)
  })
  
  return(result)
})

da = do.call(rbind, da_results_list)
saveRDS(da, "results_DA_DESeq2.rds")
write.table(da, "results_DA_DESeq2.tsv", sep='\t', quote=F, row.names=F, col.names=T)

message("DONE")

########################################################

plot_bubble <- function(df, cell_counts_df, filename) {
    plot_sig_data <- df %>%
        filter(p_val_adj < 0.05) %>%
        group_by(CellType, cluster) %>%
        summarise(n_sig_peaks = n(), .groups = 'drop')

    plot_data <- plot_sig_data %>%
        left_join(cell_counts_df, by = c("CellType", "cluster")) %>%
        mutate(n_sig_peaks = replace_na(n_sig_peaks, 0))

    p <- ggplot(plot_data, aes(x = CellType, y = cluster)) +
        geom_point(aes(size=n_cells, color=n_sig_peaks), alpha = 0.8) +
        scale_size_continuous(range = c(2, 15), name = "Number of Cells") +
        scale_color_viridis_c(name = "Number of DA Peaks\n(p.adj < 0.05)", option = "plasma") +
        labs(title = "Number of DA Peaks by Cell Type and Perturbation", x = "Cell Type", y = "Assigned perturbation") +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", hjust = 0.5))

    pdf(filename, width = 8, height = 6)
    print(p)
    dev.off()
}


plot_volcano <- function(df, filename, ct.nam) {
    p_val_thresh <- 0.05

    df <- df %>%
        mutate(
            Significance = case_when(
                p_val_adj < p_val_thresh & avg_log2FC > 0 ~ "Up-regulated (DA)",
                p_val_adj < p_val_thresh & avg_log2FC < 0 ~ "Down-regulated (DA)",
                TRUE ~ "Not significant"
            ),
            log10_p_val_adj = -log10(p_val_adj)
        )

    # Prepare data for labeling the top 10 most significant peaks per cluster
    label_df <- df %>%
        filter(Significance != "Not significant") %>%
        group_by(cluster) %>%
        slice_max(order_by = log10_p_val_adj, n = 30) %>%
        ungroup()

    p <- ggplot(df, aes(x = avg_log2FC, y = log10_p_val_adj)) +
        geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = -log10(p_val_thresh), linetype = "dashed", color = "grey50") +
        geom_text_repel(data = label_df, aes(label = gene_name), size = 3, box.padding = 0.5, point.padding = 0.5) +
        facet_wrap(~cluster, scales = "free_x") +
        scale_color_manual(values=c("Up-regulated (DA)"="firebrick", "Down-regulated (DA)"="royalblue", "Not significant"="gray70")) +
        labs(
            title = paste(glue("{ct.nam} - DA Peaks")),
            x = expression(log[2]~"Fold Change"),
            y = expression(-log[10]~"Padj"),
            color = "Peak") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

        pdf(filename, width = 10, height = 8)
        print(p)
        dev.off()
}

