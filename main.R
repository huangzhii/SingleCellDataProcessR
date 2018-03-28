# 03/28/2018 Zhi Huang
# Cell Range R Kit https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/rkit
# installation: just one line: source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
library(cellrangerRkit)
packageVersion("cellrangerRkit")
workdir = "~/Desktop/SingleCellDataProcessR/Data/"
# main starts
datapath = c("NS_084_Yoder/CB-ECFC", "NS_084_Yoder/iPS-ECFC-2", "NS_085_Yang/KO-6", "NS_085_Yang/WT-5")
for(i in 1:4){
  cellranger_pipestance_path <- paste(workdir, datapath[i], sep="")
  gbm <- load_cellranger_matrix(cellranger_pipestance_path)
  # could not find /outs/metrics_summary.csv
  dim(exprs(gbm)) # expression matrix
  dim(fData(gbm)) # data frame of genes
  dim(pData(gbm)) # data frame of cell barcodes
  expr_full_mat = as.matrix(exprs(gbm))
  write.table(expr_full_mat, file = paste(workdir, datapath[i],"_expr.csv", sep=""),row.names=T, na="",col.names=NA, sep=",")
  write.table(fData(gbm), file = paste(workdir, datapath[i],"_genes.csv", sep=""),row.names=T, na="",col.names=NA, sep=",")
  write.table(pData(gbm), file = paste(workdir, datapath[i],"_barcodes.csv", sep=""),row.names=T, na="",col.names=NA, sep=",")
}

# analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)
# 
# # t-SNE
# tsne_proj <- analysis_results$tsne
# visualize_umi_counts(gbm,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(3,4),marker_size=0.05)
