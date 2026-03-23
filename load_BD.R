library(Matrix)
source("~/kwanho/src/seurat_tools.R")
options(Seurat.object.assay.version = "v3")



# Reads in CellLevelQC output from a path
# Appends to the Seurat metadata
add_cell_level_qc <- function(qc.dir, seur.obj, base='total') {
        qc.files = list.files(qc.dir, pattern="^qc_")
	nams = gsub("^qc_|_RNA.txt$", "", qc.files)
        out = list()
        for (i in 1:length(qc.files)) {
                qc.data = read.table(file.path(qc.dir,qc.files[i]), sep='\t', row.names=1, header=T)
                rownames(qc.data) = paste0(nams[i], '_', rownames(qc.data))
                use.cbc = intersect(colnames(seur.obj), rownames(qc.data))
                qc.data = qc.data[use.cbc,]
                rm.cols = qc.data %>% summarise_all(var) %>% select_if(function(.) . == 0) %>% names()  # remove 0 variance columns
                qc.data = qc.data[,setdiff(colnames(qc.data), rm.cols)]
                for (col.nam in setdiff(colnames(qc.data), base)) {
                        qc.data[[col.nam]] = qc.data[[col.nam]]/qc.data[[base]]
                }
                out[[i]] = qc.data
        }
        qc = do.call(rbind, out)
        qc = qc[colnames(seur.obj),]
        seur.obj <- AddMetaData(seur.obj, qc)
        return(seur.obj)
}


# Function to draw scatter plot
my_qc_scatter <- function(dat, x_col, y_col, color.by, shape.by, alpha, shapes, pt.size) {
	p = ggplot(dat, aes_string(x_col, y_col, color=color.by, shape=shape.by)) +
        	geom_point(size=pt.size, alpha=alpha) +
        	scale_shape_manual(values=shapes) +
        	theme_minimal() +
        	theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
	return(p)
}


# Make plots for all cell level QC
make_qc_plots <- function(seur, group.by, shape.by, out.suffix, alpha=0.5, shapes=c(4, 19), pt.size=1,violin.width=5,add.hybrid_score=F) {
	qc.cell.lv = c('antisense','intergenic','intronic','exonic','multi','unmapped','spliced')
	xs = c(rep('nCount_RNA', 3), rep("total", length(qc.cell.lv)))
	ys = c('nFeature_RNA','percent.mito','percent.ribo',qc.cell.lv)
	if (add.hybrid_score) {
		xs = c(rep('nCount_RNA', 4), rep("total", length(qc.cell.lv)))
		ys = c('nFeature_RNA','percent.mito','percent.ribo','hybrid_score',qc.cell.lv)
	}
	
	pl = list()
	for (i in 1:length(xs)) {
        	pl[[i]] = my_qc_scatter(seur@meta.data, xs[i], ys[i], group.by, shape.by, alpha, shapes, pt.size)
	        if (i < length(xs)) {pl[[i]] = pl[[i]] + NoLegend()}
	}
	png(paste0("qc_scatter_", out.suffix, ".png"), height=12, width=10, res=300, units='in')
	print(wrap_plots(pl, ncol=3))
	dev.off()

	qc.cols = c('nCount_RNA', ys, 'total')
	MyViolinPlot(seur, features=qc.cols, group.by=group.by,
        	wid=violin.width, hei=35, ncol=1, filename=paste0("qc_violin_", out.suffix, ".png"))
}


# bd.dir = path to a directory that contains BD Pipeline output(s).
#          Each output must contain "unfiltered_RNA" directory with unzipped contents of RSEC_MolsPerCell_Unfiltered_MEX.zip
# qc.dir = path to CellLevel QC output directory. Sample names must match. Looks for "qc_" prefix.
# bd.list = alternative to 'bd.dir', can accept list of paths
LoadBD <- function(bd.dir, qc.dir, proj.nam, grnas, bd.list=c(), save.init=T, init.seur.file="seur_init.qs", 
	min.cells=10, min.features=500, min.counts=1000, npcs=30, clust.res=0.5, out.seur.file="seur_init_clustered.qs",
	add.cell.level.qc=T, add.BD.pred.cells=T, add.sample.tag=F)
{
	require(qs)
	require(Seurat)

	if (endsWith(bd.dir, "/")) bd.dir = sub("/$", "", bd.dir)
	if (endsWith(qc.dir, "/")) qc.dir = sub("/$", "", qc.dir)

	# get list of BD output directories
	print("Fetching files!")
        if(length(bd.list)==0) {
                print(paste("ls ",bd.dir,"/*/unfiltered_RNA | grep : | sed 's/://g'",sep=""))
                bd.list=system(paste("ls ",bd.dir,"/*/unfiltered_RNA | grep : | sed 's/://g'",sep=""),intern=T)
        }
	samp.nams = gsub(paste0(bd.dir,'/|/unfiltered_RNA'), "", bd.list)
	names(bd.list) = samp.nams
	print(bd.list)

	# Load counts matrix
	print("Read in!")
	seur.list <- lapply(1:length(bd.list), function(i) {
		dir = bd.list[[i]]
		nam = names(bd.list)[i]
		print(nam)
		mat <- Read10X(dir)
		print(dim(mat))
		se = CreateSeuratObject(counts=mat, project=proj.nam, min.cells=min.cells, min.features=200)
		se$sample = nam
		return(se)
	})
	names(seur.list) = names(bd.list)
	
	# Create a single Seurat object
	if (length(seur.list)>1) {
		seur <- merge(x=seur.list[[1]], y=seur.list[-1], add.cell.ids=samp.nams)
	} else {
		seur <- seur.list[[1]]
		colnames(seur) <- paste0(names(seur.list), '_', colnames(seur))
	}

	# separate gRNA and mRNA in the counts matrix
	meta = seur@meta.data
	mat = seur@assays$RNA@counts
	grna.mat = mat[grnas,]
	new.mat = mat[-which(rownames(mat) %in% grnas),]
	seur <- CreateSeuratObject(counts=new.mat, meta.data=meta, min.cells=min.cells, min.features=min.features)
	grna.mat = grna.mat[,colnames(seur)]
	seur[['gRNA']] = CreateAssayObject(counts=grna.mat, min.cells=min.cells, min.features=0)
	print(seur)

	# save init
	if (save.init) qsave(seur, init.seur.file)

	# Initial QC filtering
	seur <- subset(seur, subset = nFeature_RNA > min.features & nCount_RNA > min.counts)

	### prepare metadata
	# BD predicted cells
	if (add.BD.pred.cells) {
		bd.cells.files = gsub("unfiltered_RNA", "filtered_RNA/barcodes.tsv.gz", bd.list)
		bd.pred.cells = c()
		for (i in 1:length(bd.cells.files)) {
			f = bd.cells.files[[i]]
			nam = names(bd.cells.files)[i]
			cur.cells = paste0(nam, '_', as.character(unname(unlist(read.csv(f)))))
			bd.pred.cells = c(bd.pred.cells, cur.cells)
		}
		seur$is.cell.BD = F
		seur$is.cell.BD[bd.pred.cells] = T
	}

	if (add.sample.tag) {
		sample.tag.files = paste0(bd.dir, "/", samp.nams, "/", sub("_", "", samp.nams), "_RNA_Sample_Tag_Calls.csv")
		names(sample.tag.files) = samp.nams
		sample.tags.list = list()
		for (i in 1:length(sample.tag.files)) {
			f = sample.tag.files[[i]]
			nam = names(sample.tag.files)[i]
			cur.tags = read.csv(f, comment.char='#')
			cur.tags$Cell_Index = paste0(nam, '_', cur.tags$Cell_Index)
			sample.tags.list[[nam]] = cur.tags
		}
		sample.tags = do.call(rbind, sample.tags.list)
		seur$sample.tag = NA
		seur$sample.tag[sample.tags$Cell_Index] = sample.tags$Sample_Tag
		seur$sample.tag = paste0(seur$sample, '_', seur$sample.tag)
	}

	# Add cell level QC
	if (add.cell.level.qc) {
		seur <- add_cell_level_qc(qc.dir, seur, base='total')
	}

	# Add other QC
	seur[['percent.mito']] = PercentageFeatureSet(seur, pattern = "^mt-")
	seur[['percent.ribo']] = PercentageFeatureSet(seur, pattern = "^Rp[sl]")

	# Initial QC plots
	Idents(seur) = 'sample'
	if (add.BD.pred.cells) {
		make_qc_plots(seur, group.by='sample', shape.by='is.cell.BD', out.suffix="init")
	} else {
		make_qc_plots(seur, group.by='sample', out.suffix="init")
	}
	
	# Process data
	seur = NormalizeData(seur)
	seur<-FindVariableFeatures(seur)
	seur<-ScaleData(seur,vars.to.regress=c("nCount_RNA","percent.mito"))
	seur<-RunPCA(seur,npcs=npcs,verbose=F)
	
	# PCA plots
	pdf('pca_plots.pdf')
	print(VizDimLoadings(object = seur, dims = 1:2))
	print(DimPlot(object = seur))
	DimHeatmap(object = seur, dims = 1:9, cells = 500, balanced = TRUE)
	print(ElbowPlot(object = seur, npcs))
	dev.off()
	
	# Dimensionality reduction and clustering
	seur <- RunUMAP(seur, dims=1:npcs)
	seur <- FindNeighbors(seur, dims=1:npcs)
	seur <- FindClusters(seur, resolution=clust.res)

	# save
	qsave(seur, out.seur.file)

	# plots
	plot.feas=c('nCount_RNA','nFeature_RNA','percent.mito','percent.ribo',
		    'antisense','intergenic',"intronic","exonic","unmapped","spliced","total")
	dim.feas = c('sample','seurat_clusters','is.cell.BD')
	if (!add.cell.level.qc) {
		plot.feas=c('nCount_RNA','nFeature_RNA','percent.mito','percent.ribo')
		dim.feas = c('sample','seurat_clusters')
	}
	MyDimPlotMulti(seur, features=dim.feas, filename="umap_init_clusters.png", size=8)
	plot_feature2(seur, features=plot.feas, size=5, alpha=.5, res=500, filename="qc_feature_init_clusters.png")	
	MyViolinPlot(seur, features=plot.feas, group.by='seurat_clusters',
                wid=9, hei=36, ncol=1, filename="qc_violin_init_clusters.png")

}




