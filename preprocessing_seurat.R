# 05/04/2018 Zhi Huang
library(Seurat)
library(dplyr)
workdir = "~/Documents/DrZhiHanResearch/singlecelldata/results/new results (10X_scRNA)/"
setwd("~/Documents/DrZhiHanResearch/SingleCellDataProcessR/")

# Load the PBMC dataset
#CB-ECFC
pbmc.data <- Read10X(data.dir = "~/Documents/DrZhiHanResearch/singlecelldata/10X_scRNA/NS_084_Yoder/CB-ECFC/outs/filtered_gene_bc_matrices/GRCh38")
female = read.table("~/Documents/DrZhiHanResearch/singlecelldata/10X_scRNA/NS_084_Yoder/NS_084_Yoder_CB-ECFC_genderClassfication_female.txt", header = T)
#iPS-ECFC-2
pbmc.data <- Read10X(data.dir = "~/Documents/DrZhiHanResearch/singlecelldata/10X_scRNA/NS_084_Yoder/iPS-ECFC-2/outs/filtered_gene_bc_matrices/GRCh38")
female = read.table("~/Documents/DrZhiHanResearch/singlecelldata/10X_scRNA/NS_084_Yoder/NS_084_Yoder_iPS-ECFC_genderClassfication_female.txt", header = T)
#KO-6
pbmc.data <- Read10X(data.dir = "~/Documents/DrZhiHanResearch/singlecelldata/10X_scRNA/NS_085_Yang/KO-6/outs/filtered_gene_bc_matrices/GRCh38")
#WT-5
pbmc.data <- Read10X(data.dir = "~/Documents/DrZhiHanResearch/singlecelldata/10X_scRNA/NS_085_Yang/WT-5/outs/filtered_gene_bc_matrices/GRCh38")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size
sparse.size <- object.size(x = pbmc.data)
sparse.size
dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 5, min.genes = 200, 
                           project = "10X_PBMC")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(3000, -Inf), high.thresholds = c(9500, 0.1),
                    cells.use = female[,1])

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)


# Scaling the data and removing unwanted sources of variation
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

# Perform linear dimensional reduction
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
