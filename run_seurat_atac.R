# conda activate multiome

library(Signac)
library(Seurat)
library(dplyr)
library(GenomicRanges)
library(AnnotationHub)
library(ggplot2)
library(patchwork)
library(qs)
library(glue)
source("~/kwanho/src/density_scatter.R")
source("~/kwanho/src/seurat_tools.R")
options(Seurat.object.assay.version = "v3")


############################################################

outdir = "/stanley/levin_xinjin/kwanho/BD_Rhapsody/analysis/260121_in_vivo/BD_Pipeline"
pattern = "^vivo"

bd.samps = list.files(outdir, pattern=pattern)
bd.dirs = file.path(outdir, bd.samps)

############################################################

print("Start!")
print(glue("BD pipeline path: {outdir}"))
print(glue("Sample: {bd.samps}"))

# Read in
obj.list <- lapply(1:length(bd.dirs), function(i) {
	cur.dir = bd.dirs[i]
	cur.samp = bd.samps[i]
	print(cur.dir)
	cts <- ReadMtx(mtx=glue("{cur.dir}/filtered_ATAC/atac-matrix.mtx.gz"),
                 cells=glue("{cur.dir}/filtered_ATAC/atac-barcodes.tsv.gz"),
                 features=glue("{cur.dir}/filtered_ATAC/atac-features.tsv.gz"))
	frag.path <- glue("{cur.dir}/{sub('_','',cur.samp)}_RNA_ATAC_Fragments.bed.gz")
	chrom.assay <- CreateChromatinAssay(counts=cts, sep = c(":", "-"), fragments = frag.path, genome = 'mm39')
	cur.obj = CreateSeuratObject(counts = chrom.assay, project='Xinhe_in_vivo_ATAC', assay='ATAC')
	return(cur.obj)
})

qsave(obj.list, "obj_list_init.qs")

# common peak set
combined.peaks <- UnifyPeaks(object.list=obj.list, mode="reduce")

# cells
wl.files = list.files(pattern="^bc_whitelist")

# common peak counts
peak.cts.list <- lapply(1:length(wl.files), function(i) {
	wl.file = wl.files[i]
	obj = obj.list[[i]]
	frag=obj@assays$ATAC@fragments[[1]]
	FeatureMatrix(fragments=frag, features=combined.peaks, sep=c(":", "-"), cells=scan(wl.file))
})
qsave(peak.cts.list, "common_peak_counts_list.qs")

# merge peak counts
all.peak.cts <- do.call(cbind, peak.cts.list)
colnames(all.peak.cts) <- unlist(lapply(1:length(bd.samps), function(i) {
	paste0(bd.samps[i], "_", colnames(peak.cts.list[[i]]))
}))
qsave(all.peak.cts, "peak_counts_combined.qs")

# create Seurat object
chrom_assay = CreateChromatinAssay(counts=all.peak.cts, sep=c(":", "-"), genome="mm39")
seur <- CreateSeuratObject(counts=chrom_assay, assay="peaks", project="Xinhe_in_vivo_ATAC")
DefaultAssay(seur) <- "peaks"

# save
qsave(seur, "seur_atac_init.qs")

###################################################################################
if (F) {

# Add gene annotation
hub <- AnnotationHub()
query_results <- query(hub, c("EnsDb", "Mus musculus", "GRCm39"))
print(tail(query_results))
latest_id <- "AH119358"
ensdb = hub[[latest_id]]
# Extract gene, transcript, and exon features separately
annotations_genes <- genes(ensdb)
annotations_transcripts <- transcripts(ensdb)
annotations_exons <- exons(ensdb)
# Add metadata to each feature type before combining them
mcols(annotations_genes)$type <- "gene"
mcols(annotations_genes)$gene_name <- annotations_genes$gene_name
mcols(annotations_genes)$gene_biotype <- annotations_genes$gene_biotype
# Match transcripts to genes
transcript_info <- ensembldb::select(
  ensdb,
  keys = annotations_transcripts$tx_id,
  keytype = "TXID",
  columns = c("GENEID", "GENENAME", "GENEBIOTYPE")
)
m <- match(annotations_transcripts$tx_id, transcript_info$TXID)
mcols(annotations_transcripts)$type <- "transcript"
mcols(annotations_transcripts)$gene_id <- transcript_info$GENEID[m]
mcols(annotations_transcripts)$gene_name <- transcript_info$GENENAME[m]
mcols(annotations_transcripts)$gene_biotype <- transcript_info$GENEBIOTYPE[m]
# Match exons to genes
exon_info <- ensembldb::select(
  ensdb,
  keys = annotations_exons$exon_id,
  keytype = "EXONID",
  columns = c("GENEID", "GENENAME", "GENEBIOTYPE")
)
m <- match(annotations_exons$exon_id, exon_info$EXONID)
mcols(annotations_exons)$type <- "exon"
mcols(annotations_exons)$gene_id <- exon_info$GENEID[m]
mcols(annotations_exons)$gene_name <- exon_info$GENENAME[m]
mcols(annotations_exons)$gene_biotype <- exon_info$GENEBIOTYPE[m]
# combine gene-/transcript-/exon-level annotations
annotations <- c(annotations_genes, annotations_transcripts, annotations_exons)
annotations <- keepStandardChromosomes(annotations, species = "Mus musculus", pruning.mode = "coarse")
seqlevelsStyle(annotations) <- "UCSC"  # 'chr' prefix in chromsome names to match those in fragments file
Annotation(seur) <- annotations

###################################################################################

# Add combined fragments file
frag.file = "../BD_Pipeline/fragments.bed.gz"
frags <- CreateFragmentObject(path=frag.file)
Fragments(seur[["peaks"]]) <- frags

qsave(seur, "seur_atac.qs")

# QC
seur <- NucleosomeSignal(object = seur)  # add nucleosome signal
seur <- TSSEnrichment(object = seur)  # add TSS enrichment
all_counts <- CountFragments(fragments=frag.file)
rownames(all_counts) <- all_counts$CB
all_counts$total_fragments <- all_counts[,"frequency_count"]
all_counts = all_counts[colnames(seur),]
seur <- AddMetaData(seur, all_counts)
seur <- FRiP(seur, assay='peaks', total.fragments="total_fragments", col.name="FRiP_reads_in_peaks")  # add FRiP
blacklist_mm39 = readRDS("mm39.excluderanges.rds")
seur$blacklist_ratio <- FractionCountsInRegion(object = seur, assay = 'peaks',  regions = blacklist_mm39)  # add blacklist ratio

# save
qsave(seur, "seur_atac_qc.qs")

# QC plots
seur$sample = apply(stringr::str_split(colnames(seur), '_', simplify=T)[,1:2], 1, paste, collapse='_')
MyViolinPlot(seur, features=c('nCount_peaks','TSS.enrichment','blacklist_ratio','nucleosome_signal','FRiP_reads_in_peaks'),
	filename="qc_atac_violin.pdf", nc=1, hei=20, wid=12, add.median=T, x.lab="", group.by='sample')
p1=ggplot(seur@meta.data, aes(nCount_peaks, TSS.enrichment)) + geom_point() + scale_x_log10() +
	geom_hline(yintercept=1) + geom_vline(xintercept=200)
p2=ggplot(seur@meta.data, aes(nCount_peaks, nucleosome_signal)) + geom_point() + scale_x_log10() +
	geom_hline(yintercept=0.3) + geom_vline(xintercept=200)
p3=ggplot(seur@meta.data, aes(nCount_peaks, blacklist_ratio)) + geom_point() + scale_x_log10() +
	geom_hline(yintercept=0.075) + geom_vline(xintercept=200)
p4=ggplot(seur@meta.data, aes(nCount_peaks, FRiP_reads_in_peaks)) + geom_point() + scale_x_log10() +
	geom_hline(yintercept=0.1) + geom_vline(xintercept=200)
p=wrap_plots(list(p1,p2,p3,p4), ncol=1)
ggsave('qc_atac_scatter_qc.png', plot=p, height=15, width=7)
DensityScatter(seur, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave("qc_atac_density_scatter_nPeaks_TSSenrich.png", height=5, width=7)
seur$nucleosome_group <- ifelse(seur$nucleosome_signal>0.2, 'NS>0.2', 'NS<0.2')
seur$nucleosome_group <- paste0(seur$sample, '_', seur$nucleosome_group)
FragmentHistogram(object = seur, region = "chr1-3000000-5000000", group.by='nucleosome_group')
ggsave("qc_hist_frag_length_NS_signal.png", height=5, width=7)
FragmentHistogram(object = seur, region = "chr1-3000000-5000000", group.by='sample')
ggsave("qc_hist_frag_length_samples.png", height=5, width=10)

}
