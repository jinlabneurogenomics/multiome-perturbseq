library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(glmGamPoi)
library(ggplot2)
library(cowplot)
library(patchwork)
library(qs)
library(glue)
source("~/kwanho/src/seurat_tools.R")


rna = qread("../RNA/seur_final_v2.qs")
atac = qread("../ATAC/seur_atac_qc.qs")

grnas = gsub("_", "-", scan("grna_list.txt", ''))
all_guides = c(grnas, 'Doublet','Not assigned')
n_grna = length(all_guides)-2

guide.cols = viridis::viridis(n_grna)
my.cols = c(guide.cols,'gray80','gray95')
names(my.cols) = all_guides

my.cols2 = c(Singlet="black", Doublet="gray70", "Not assigned"="gray95")

ct.cols <- c(
  "Car3"    = "#E38B5F",
  "L2-3.IT" = "#2F7EB8",
  "L4-5.IT" = "#9DCEE5",
  "L5.IT"   = "#5E9FB3",
  "L5.NP"   = "#B88C6A",
  "L5.ET"   = "#E2C45D",
  "L6.IT"   = "#76B7D8",
  "L6.CT"   = "#B8D09A",
  "L6b"     = "#6E9E86",
  "Pvalb"   = "#B6A0D6",
  "Sst"     = "#D3A5C4",
  "Vip"     = "#F59AAE"
)

##########################################

# subset to common cells
common.cells = intersect(colnames(rna), colnames(atac))  # 84887 cells
#rna = subset(rna, cells=common.cells)
atac = subset(atac, cells=common.cells)

# integrate data
seur <- rna
seur[['peaks']] = atac@assays$peaks
atac.meta = atac@meta.data[,c('nCount_peaks','nFeature_peaks',
			      'nucleosome_signal','nucleosome_percentile','TSS.enrichment','TSS.percentile','total_fragments',
			      'mononucleosomal','nucleosome_free','reads_count','FRiP_reads_in_peaks','blacklist_ratio')]
seur <- AddMetaData(seur, atac.meta)

qsave(seur, "seur_multiome_init.qs")

# RNA analysis
#DefaultAssay(seur) <- "RNA"
#seur <- SCTransform(seur, vars.to.regress=c('percent.mito','intronic'), verbose=FALSE) %>% 
#	RunPCA() %>% 
#	RunUMAP(dims=1:30, reduction.name='umap.rna', reduction.key='rnaUMAP_')
#
#qsave(seur, "seur_multiome_1.qs")

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(seur) <- "peaks"
seur <- RunTFIDF(seur)
seur <- FindTopFeatures(seur, min.cutoff='q0')
seur <- RunSVD(seur)
seur <- RunUMAP(seur, reduction='lsi', dims=2:30, reduction.name="umap.atac", reduction.key="atacUMAP_")

qsave(seur, "seur_multiome_2.qs")

#DepthCor(seur)
#ggsave("qc_depth_corr.png", height=5, width=7)

# WNN
seur <- FindMultiModalNeighbors(seur, reduction.list=list("harmony", "lsi"), dims.list=list(1:50, 2:50))
seur <- RunUMAP(seur, nn.name="weighted.nn", reduction.name="wnn.umap", reduction.key="wnnUMAP_")
seur <- FindClusters(seur, graph.name="wsnn", algorithm=3, verbose=FALSE, resolution=.3)

levels(seur$CellType) = gsub("_", ".", levels(seur$CellType))

qsave(seur, "seur_multiome_final_harmony.qs")

###################################################################

# visualize
p1 <- DimPlot(seur, reduction="umap", group.by="CellType", label=TRUE, label.size=3, repel=TRUE, split.by='sample', raster=T,
	raster.dpi=c(1024,1024), pt.size=3, cols=ct.cols) +
	ggtitle("RNA") + NoLegend()
p2 <- DimPlot(seur, reduction="umap.atac", group.by="CellType", label=TRUE, label.size=3, repel=TRUE, split.by='sample', raster=T,
        raster.dpi=c(1024,1024), pt.size=3, cols=ct.cols) +
	ggtitle("ATAC") + NoLegend()
p3 <- DimPlot(seur, reduction="wnn.umap", group.by="CellType", label=TRUE, label.size=3, repel=TRUE, split.by='sample', raster=T,
        raster.dpi=c(1024,1024), pt.size=3, cols=ct.cols) +
	ggtitle("WNN")
p4 <- DimPlot(seur, reduction="wnn.umap", group.by="assign", label=TRUE, label.size=3, repel=TRUE, split.by='sample', raster=T,
        raster.dpi=c(1024,1024), pt.size=3, cols=my.cols) +
	ggtitle("WNN - gRNA assignment")
p5 <- DimPlot(seur, reduction="wnn.umap", group.by="assign_broad", label=F, split.by='sample', cols=my.cols2, raster=T,
        raster.dpi=c(1024,1024), pt.size=3) +
	ggtitle("WNN - gRNA assignment (broad)")
wrap_plots(list(p1, p2, p3, p4, p5), ncol=1)
ggsave("umap_multiome_harmony.pdf", width=25, height=25)

DefaultAssay(seur) = 'RNA'
plot.feas = c('nCount_RNA','nCount_peaks','nFeature_RNA','nFeature_peaks','antisense','intergenic','intronic',
                'hybrid_score','percent.mito','percent.ribo','nucleosome_signal','TSS.enrichment','mononucleosomal',
                'nucleosome_free','FRiP_reads_in_peaks','blacklist_ratio')
plot_feature2(seur, features=plot.feas, nc=4, filename="qc_feature_plot_final_harmony.png", size=4, plot.width=30, 
	plot.height=15, alpha=.5, reduction='wnn.umap')
MyViolinPlot(seur, features=plot.feas, filename="qc_violin_plot_CellType_final.png", group.by='CellType', rm.legend=F, 
	add.median=T, height.multiplier=1.25, width.multiplier=1.5, median.pos=position_dodge(width=1))
MyViolinPlot(seur, features=plot.feas, filename="qc_violin_plot_samples_final.png", group.by='sample', rm.legend=F, 
	add.median=T, height.multiplier=1.25, width.multiplier=1.5, median.pos=position_dodge(width=1))


