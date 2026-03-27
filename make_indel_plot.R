library(Seurat)
library(dplyr)
library(qs)
library(ggplot2)
library(scales)
library(patchwork)
library(stringr)
library(tidyr)
library(glue)
library(viridis)


create_indel_bubble_plot <- function(data, clust_order=NULL, guide_order=NULL, min_reads=0, hline_y=8,
                                            cols=rev(magma(256)[2:220]), scale_rows=TRUE) {
  
  plot_data <- data %>%
    filter(!is.na(Reads) & Reads > min_reads) %>%
    mutate(
      INDEL_Ratio = replace_na(INDEL_Ratio, 0),
      Guide = if (!is.null(guide_order)) factor(Guide, levels = guide_order) else as.factor(Guide),
      Cluster = if (!is.null(clust_order)) factor(Cluster, levels = clust_order) else as.factor(Cluster)
    )

  if (scale_rows) {
    plot_data <- plot_data %>%
      group_by(Guide) %>%
      mutate(
        INDEL_Ratio = (INDEL_Ratio - min(INDEL_Ratio, na.rm = TRUE)) / (max(INDEL_Ratio, na.rm = TRUE) - min(INDEL_Ratio, na.rm = TRUE))
      ) %>%
      ungroup()
  }

  if (nrow(plot_data) == 0) {
    message("No data points with Reads > 0 to plot.")
    return(invisible(NULL))
  }
  
  # Determine the title and legend name based on the scaling
  title_text <- ifelse(scale_rows, "Scaled INDEL Ratio and Read Count by Guide and Cluster", 
	"INDEL Ratio and Read Count by Guide and Cluster")
  color_name <- ifelse(scale_rows, "Scaled INDEL Ratio", "INDEL Ratio")

  p <- ggplot(plot_data, aes(x = Cluster, y = Guide, color = INDEL_Ratio, size = Reads)) +
    geom_point(alpha = 0.8) +
    #geom_hline(yintercept=hline_y+0.5, linetype = "dashed", color = "black") +
    #facet_wrap(~ Method, ncol = 2) +
    scale_color_gradientn(name = color_name, colors=cols, limits=c(0, NA)) +
    scale_size_continuous(
      name = "Total Reads",
      range = c(1, 20),
      breaks = c(1, 10, 100, 1000, 5000),
      limits = c(1, 5000)
    ) +
    labs(
      title = title_text,
      subtitle = glue("Min total reads = {min_reads}"),
      x = "Cells assigned to a gRNA",
      y = "gRNA"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12, face = "bold", margin = margin(t = 5, b = 5)),
      strip.background = element_rect(fill = "#f0f0f0", color = "grey50"),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = margin(t = 1, r = 2, b = 1, l = 2, unit = "cm")
    )

  return(p)
}

#res <- read.csv("results_full.csv")
res <- read.csv("results_fully_spanning.csv")
res$indel = res$Insert + res$Delete
res$sample = str_extract(res$File, pattern="vivo[1-5]")
res$sample = gsub("vivo", "vivo_", res$sample)
res$Cell = paste0(res$sample, '_', res$CBC)

seur <- qread("../../RNA/seur_final_v2.qs")

cluster_order <- scan("../../Seurat/grna_list.txt","")
guides = gsub("-", "_", cluster_order)

# high coverage region targeting guides
#guides_hc=c('494_Ptprd_g1A','494_Ptprd_g2B','494_Rtn1_g2A','494_Rtn1_g2B','496_Celf2_g2A','496_Celf2_g2B','496_Nrxn1_g1A','496_Nrxn1_g2B')
#guides = c(guides_hc, setdiff(names(table(res$Guide)), guides_hc))

for (ch in 1:5) {
	print(ch)
	out.li = list()
	met = paste0("vivo_", ch)
        sres <- res %>% filter(grepl(met, sample))
        for (guide in guides) {
                print(guide)
                li = list()
                dat = sres %>% filter(Guide == guide)
                print(dim(dat))
                if (nrow(dat)>0) {
                        for (pert in setdiff(names(table(seur[['assign']])), c('Doublet','Not assigned'))) {
                                cur_assign = seur[['assign']]
                                cells = rownames(cur_assign)[grepl(pert, cur_assign[,1])]
                                cur_dat = dat %>% filter(Cell %in% cells)
                                li[[paste("assign", pert)]] = c(colSums(cur_dat[,2:4]), "Cluster"=pert)
                        }
                	out = as.data.frame(do.call(rbind, li))
                	cols.num = c('Insert','Delete','Reads')
                	out[cols.num] <- sapply(out[cols.num], as.numeric)
	                out = out %>% as.data.frame() %>% mutate(INDEL = Insert + Delete) %>%
                        	mutate(INDEL_Ratio=INDEL/Reads, Guide=guide)
                	print(out)
                	out.li[[guide]] = out %>% as_tibble()
                	#write.table(out, paste0("results_INDEL_", guide, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
		}
        }
	pdat = do.call(rbind, out.li)
	write.table(pdat, glue("table_INDEL_vivo{ch}.tsv"), sep='\t', quote=F, row.names=F, col.names=T)
	create_indel_bubble_plot(pdat, cluster_order, guides, min_reads=5)
	#create_indel_bubble_plot(pdat, cluster_order, guides_hc, min_reads=5)
	ggsave(glue("bubble_INDEL_ratio_vivo{ch}_min5.pdf"), width=10, height=13)
}


