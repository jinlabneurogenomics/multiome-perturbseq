library(Matrix)
library(corrplot)
library(tidyr)
library(scales)
source("~/kwanho/src/seurat_tools.R")
source("~/kwanho/src/load_BD.R")
options(Seurat.object.assay.version = "v3")


############################################################
# Helper functions
# Function to compute correlation of gRNAs
# cts = matrix of gRNAs X CBC
get_corr <- function(cts, outnam="") {
	dat = t(cts)
	cor.val = cor(as.matrix(dat), method='spearman')
	pdf(paste0("guide_correlation_", outnam, ".pdf"))
	print(corrplot(cor.val, tl.col='black',title=paste("BD Spearman correlation", outnam),mar=c(0,0,4,0), 
		tl.offset=1, col=rev(COL2("RdBu", 200))))
	dev.off()
}

# Subset seurat object, reprocess, and make plots
subset_reprocess <- function(seur.obj, rm.cls, out.suffix, dim.plot.feas, plot.feas, group.by, shape.by) {
	if (length(rm.cls)>0) {
		print("subset!")
		seur.obj <- subset(seur.obj, idents=rm.cls, invert=T)
	}
	print(seur.obj)
	seur.obj <- StandardSeuratProcessing(seur.obj, filename=glue("seur_{out.suffix}.qs"), is.mouse=T, 
					 out.name=out.suffix, filter.cells=F, save.plots=F)

	# Plots
	MyDimPlotMulti(seur.obj, features=dim.feas, filename=glue("umap_{out.suffix}.png"), size=8)
	plot_feature2(seur.obj, features=plot.feas, size=5, alpha=.5, res=500, filename=glue("qc_feature_{out.suffix}.png"))
	make_qc_plots(seur.obj, group.by, shape.by, out.suffix, violin.width=7)
	return(seur.obj)
}

# get second max value from counts matrix
get_second_max <- function(x) {
  sorted_values <- sort(x, decreasing = TRUE)
  if (length(sorted_values) < 2) {
    return(0) # Handle cases with fewer than 2 values
  } else {
    return(sorted_values[2])
  }
}

# Scatter plot of the top two gRNA counts in each cell
scatter_top2_guides <- function(x, filename, assign.col, order=NULL, merge_singlet=T, not.singlet=c('Doublet','Not assigned')) {
  g1 = apply(x, 1, max)
  g2 = apply(x, 1, get_second_max)
  df = data.frame(max=g1, second_max=g2)
  df$assign = seur[[assign.col]][rownames(x),]
  if (!is.null(order)) { df$assign = factor(df$assign, levels=order) }
  if (merge_singlet) { levels(df$assign)[which(!levels(df$assign) %in% not.singlet)] = 'Singlet'}
  ggplot(df, aes(x=max, y=second_max, color=assign)) + geom_point(size=.25, alpha=.5) + 
	geom_abline(slope=1, intercept=0, color='black',
	linetype=2) + theme_minimal()
  ggsave(filename, dpi=600, height=4, width=6)
}

############################################################

# Set up!
proj.nam = "Xinhe_in_vivo_perturb"

grnas = scan("../Seurat/grna_list.txt",'')
constant.seqs = scan("constant_seqs.txt",'')

dim.feas = c('sample','sample.tag','seurat_clusters','CellType')
plot.feas=c('nCount_RNA','nFeature_RNA','percent.mito','percent.ribo','hybrid_score',
        'antisense','intergenic',"intronic","exonic","unmapped","spliced","total")

############################################################
# combine vivo 1
vivo25 <- qread("../Seurat/seur_final.qs")
vivo1 <- qread("../Seurat/vivo1/seur_vivo1_ready.qs")

seur <- merge(vivo1, vivo25)

############################################################

print("Start!")

meta = seur@meta.data
mat = seur@assays$RNA@counts
gmat = seur@assays$gRNA@counts
new.mat = mat[-which(rownames(mat) %in% constant.seqs),]
seur <- CreateSeuratObject(counts=new.mat, meta.data=meta)
seur[['gRNA']] = CreateAssayObject(counts=gmat)

qsave(seur, "seur_final_v2_init.qs")

# Harmony
library(harmony)
seur <- seur %>% 
	FindVariableFeatures() %>% 
	ScaleData(vars.to.regress=c("nCount_RNA","percent.mito")) %>% 
	RunPCA() %>% 
	RunHarmony("sample", plot_convergence=F) %>%
	RunUMAP(reduction = "harmony", dims=1:20) %>%
	FindNeighbors(reduction="harmony", dims=1:20) %>%
	FindClusters(resolution=0.5)

qsave(seur, "seur_final_v2.qs")

# Plots
MyDimPlotMulti(seur, features=dim.feas, filename="umap_final_v2.png", size=8, show.label=c(F,F,T,T))
plot_feature2(seur, features=plot.feas, size=5, alpha=.5, res=500, filename="qc_feature_final_v2.pdf", pt.size=2,
	raster=T, raster.dpi=c(1024,1024))
make_qc_plots(seur, "CellType", "is.cell.BD", "final_v2", violin.width=7, add.hybrid_score=T)

# Check gene expression
plot.genes = c('Car3','Snap25','Cux1','Cux2','Satb2','Nectin3','Rorb','Bcl11b','Foxp2','Scnn1a','Etv1','Fezf2','Trhr','Ccn2','Parm1',
        'Gad1','Gad2','Pvalb','Sst','Vip')
plot_feature2(seur, plot.genes, size=5, alpha=.5, res=500, filename="feature_genes_final_v2.pdf")
DotPlot(seur, group.by='CellType', features=plot.genes) + scale_colour_gradient2(low = "blue", mid = "gray90", high = "red2") +
        coord_flip() + theme_light()
ggsave("dotplot_genes_final_v2.pdf", height=6, width=10, dpi=500)

############################################################
# plot samples (aggregate sample tags)

tab=as.data.frame(table(seur$sample))
colnames(tab) = c('sample','num.cells')
ggplot(tab, aes(x = sample, y = num.cells)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Sample",
       y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")
ggsave("barchart_numCells_sample.pdf")

PropsPlot(seur, my.group='CellType',my.sample='sample', name.prefix="barchart_props_sample_cluster_composition", dev='pdf')

MyViolinPlot(seur, features=plot.feas, group.by='sample',
                wid=9, hei=35, ncol=1, filename="qc_violin_vivo_sample.pdf")
MyViolinPlot(seur, features=plot.feas, group.by='CellType',
                wid=12, hei=35, ncol=1, filename="qc_violin_vivo_CellType.pdf")

############################################################
# plot gRNA assignment

all_guides = c(grnas, 'Doublet','Not assigned')
n_grna = length(all_guides)-2
mets = c('assign', 'assign_broad')
grp.ord = c('Singlet','Doublet','Not assigned')

# order guides
seur@meta.data[,mets[1]] = factor(seur@meta.data[,mets[1]], levels=all_guides)
seur@meta.data[,mets[2]] = factor(seur@meta.data[,mets[2]], levels=grp.ord)

# set colors
guide.cols = viridis::viridis(n_grna)
my.cols = c(guide.cols,'gray80','gray95')
names(my.cols) = all_guides
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

# umap
print("plot gRNA assignment in UMAP")
MyDimPlotMulti(seur, features='assign', filename="umap_guide_assign.pdf", split.by='sample', nc=1, cols=my.cols,
        n.legend.rows=5, add.plot.height=1, add.plot.width=12, show.label=F, pt.size=.1, alpha=.5, size=6)
MyDimPlotMulti(seur, features='CellType', filename="umap_cell_types.pdf", split.by='sample', nc=1, cols=ct.cols,
        n.legend.rows=5, add.plot.height=1, add.plot.width=12, show.label=F, pt.size=.1, alpha=.5, size=6)

# nDEGs
exact_counts <- out %>%
  filter(fdrtool_q < 0.05) %>%
  count(cell_type, comparison, name = "n_sig")

ggplot(exact_counts, aes(x = comparison, y = cell_type, size = n_sig)) +
  geom_point(color = "steelblue") +
  geom_text(aes(label = n_sig), size = 3, color = "white") +
  scale_size_continuous(range = c(5, 15)) +
  theme_minimal() +
  labs(size="nDEGs",
       x = "Comparison",
       y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("dotplot_nDEGs.pdf")

# num cells per pert
# broad
tab = table(seur$sample, seur$assign_broad)
write.table(tab, "table_guide_assignment_counts_broad.tsv", sep='\t', row.names=F, col.names=T, quote=F)

tab <- tab / rowSums(tab)
df <- tab %>% as.data.frame()
colnames(df) <- c('Pool','Assignment','Fraction')
ggplot(df, aes(x = Assignment, y = Fraction, fill = Assignment)) +
  geom_col() +
  facet_wrap(~Pool, ncol=3) +
  theme_minimal() +
  labs(
    x = "Assignment",
    y = "Fraction of cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("barchart_guide_assignment_broad.pdf", height=5, width=10)

# all guides
tab = table(seur$sample, seur$assign)
write.table(tab, "table_guide_assignment_counts.tsv", sep='\t', row.names=F, col.names=T, quote=F)

tab <- tab / rowSums(tab)
df <- tab %>% as.data.frame()
colnames(df) <- c('Pool','Assignment','Fraction')
ggplot(df, aes(x = Assignment, y = Fraction, fill = Assignment)) +
  geom_col() +
  facet_wrap(~Pool, ncol=3) +
  theme_minimal() +
  labs(
    x = "Assignment",
    y = "Fraction of cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("barchart_guide_assignment.pdf", height=10, width=14)

ggplot(df %>% filter(!Assignment %in% c('Doublet','Not assigned')), aes(x = Assignment, y = Fraction, fill = Assignment)) +
  geom_col() +
  facet_wrap(~Pool, ncol=3) +
  theme_minimal() +
  labs(
    x = "Assignment",
    y = "Fraction of cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("barchart_guide_assignment_singletsOnly.pdf", height=10, width=14)


# plot top two guides in cells
scatter_top2_guides(t(grna2), "plot_scatter_top2_gRNA_counts.vivo2.pdf", "assign_broad", order=grp.ord, merge_singlet=F)
scatter_top2_guides(t(grna3), "plot_scatter_top2_gRNA_counts.vivo3.pdf", "assign_broad", order=grp.ord, merge_singlet=F)
scatter_top2_guides(t(grna4), "plot_scatter_top2_gRNA_counts.vivo4.pdf", "assign_broad", order=grp.ord, merge_singlet=F)
scatter_top2_guides(t(grna5), "plot_scatter_top2_gRNA_counts.vivo5.pdf", "assign_broad", order=grp.ord, merge_singlet=F)

############################################################

# 










############################################################

# Figures


MyViolinPlot(seur, features=plot.feas, group.by='orig.ident', nc=3, width.multiplier=8, filename="violin_qc_2A_vs_2C.pdf")

seur1 = subset(seur, orig.ident=='2C')
seur2 = subset(seur, orig.ident=='2A')

myStackedVlnPlot2 <- function(seur, features, group.by, title="", pt.size=1, cols=NULL, assay='RNA', data='data', facet.scales=NULL) {
require(Seurat)
require(ggplot2)
require(cowplot)
require(patchwork)
require(dplyr)
seur@active.assay=assay
if (!all(features %in% rownames(seur))) {
        print("Following genes are not found!")
        print(features[!(features %in% rownames(seur))])
        features = features[features %in% rownames(seur)]
}
df = as.data.frame(t(as.matrix(seur@assays[[assay]]@data[features,])))
if (data=='counts') {
df = as.data.frame(t(as.matrix(seur@assays[[assay]]@counts[features,])))
}
df$Cell = rownames(df)
Idents(seur) <- group.by
df$Idents = Idents(seur)
#df$Idents = seur@meta.data[,group.by]

if (is.null(cols)) {
cols = hue_pal()(nlevels(seur))
}

dat <- reshape2::melt(df, id.vars = c("Cell","Idents"), measure.vars = features,
                variable.name = "Feat", value.name = "Expr")
dat$Idents = factor(dat$Idents, c(paste0('TE.', grnas), 'Doublet','Not assigned'))

ylims <- range(dat$Expr, na.rm = TRUE)
if (is.null(facet.scales)) facet.scales='free'

a <- ggplot(dat, aes(factor(Idents), Expr, fill = NULL)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        geom_jitter(size=pt.size, position=position_jitter(seed = 1, width = 0.2)) +
        scale_fill_manual(values=cols) +
        scale_y_continuous(expand = c(0, 0), position="right", trans = "log10") +
        facet_grid(rows = vars(Feat), scales = facet.scales, switch = "y") +
        theme_cowplot(font_size = 12) +
        #theme_classic() +
        theme(legend.position = "none", panel.spacing = unit(1, "lines"),
              #panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
        #geom_boxplot(width=0.1, outlier.shape=NA) +
        stat_summary(fun=median, geom='crossbar', color='red', width=.6, linewidth=.2) +  # width=0.1
        ggtitle(title) + xlab("Cells assigned to a gRNA") + ylab("gRNA counts")
return(a)
}

myStackedVlnPlot2(seur1, features=grnas, group.by='assign', pt.size=.1, data='counts', assay='gRNA')
ggsave("violin_gRNA_expr_2C_counts_v2.pdf", height=10, width=7)
myStackedVlnPlot2(seur2, features=grnas, group.by='assign', pt.size=.1, data='counts',assay='gRNA')
ggsave("violin_gRNA_expr_2A_counts_v2.pdf", height=10, width=7)


