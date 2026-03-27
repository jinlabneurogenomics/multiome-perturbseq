# Multiome Perturb-seq Analysis Pipeline

This repository contains the computational pipeline used to process, analyze, and visualize *in vivo* single-nucleus multiomic Perturb-seq data (BD Rhapsody). The workflow covers cell-level quality control, single-modality and joint multiomic processing, guide RNA assignment, validation of gRNA assignment with INDEL in reads, and downstream differential analysis.

## File Descriptions

* **`driver.cell_level_qc.sh`** A shell script that automates the cell-level quality control (QC) steps using [CellLevel_QC](https://github.com/seanken/CellLevel_QC).

* **`load_BD.R`** An R script with functions to load multiple outputs of BD Rhapsody Sequencing Analysis Pipeline [v3.0](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/) into a single Seurat object. The `LoadBD` function reads unfiltered data and optionally appends QC measures from CellLevel_QC.

* **`process_RNA.R`** Processes snRNA-seq data. Includes steps for RNA normalization, scaling, clustering, and cell-type annotation.

* **`process_ATAC.R`** Processes snATAC-seq data. Handles joint peak calling across samples, normalization, and dimensionality reduction.

* **`process_Multiome.R`** Integrates the RNA and ATAC modalities (Weighted Nearest Neighbor approach) to create joint multiomic embeddings.

* **`assign_guides.py`** A Python script used to assign CRISPR gRNAs to individual nuclei, establishing the perturbation identity of each cell. The assignment leveraged Gaussian mixed model-based approach implemented in [Crispat](https://github.com/velten-group/crispat).

* **`GetInsertDelete.py`** A Python script utilizing [pysam](https://pysam.readthedocs.io/en/latest/api.html) to scan BAM files for indels spanning the gRNA target regions, used to validate the guide assignment.

* **`make_indel_plot.R`** An R script that visualizes the indel frequencies and editing outcomes calculated by the `GetInsertDelete.py` script.

* **`DE_test.R`** Performs DE test (edgeR+fdrtool) on the RNA data to identify DEGs resulting from specific CRISPR perturbations.

* **`DA_test.R`** Performs DA test (edgeR+fdrtool) on the ATAC data to identify DARs resulting from specific CRISPR perturbations..
