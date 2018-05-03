# 05/03/2018 Zhi Huang
workdir = "~/Documents/DrZhiHanResearch/singlecelldata/results/new results (10X_scRNA)/"

exprs.list = c("NS_084_Yoder/CB-ECFC_expr.csv",
          "NS_084_Yoder/iPS-ECFC-2_expr.csv",
          "NS_085_Yang/KO-6_expr.csv",
          "NS_085_Yang/WT-5_expr.csv")

genes.list = c("NS_084_Yoder/CB-ECFC_genes.csv",
          "NS_084_Yoder/iPS-ECFC-2_genes.csv",
          "NS_085_Yang/KO-6_genes.csv",
          "NS_085_Yang/WT-5_genes.csv")

for (i in 1:length(genes.list)){
  
  expr = read.csv(file=paste(workdir, exprs.list[i], sep = ""), header=T, row.names = 1, sep=",")
  gene = read.csv(file=paste(workdir, genes.list[i], sep = ""), header=T, sep=",")
  
  cells = colnames(expr)
  expr.binary = expr>0
  unique.gene.counts = colSums(expr.binary)
  filter.index.1 = unique.gene.counts < 9500
  filter.index.2 = unique.gene.counts > 3000
  
  mt.gene.index = "MT-" == substr(gene[,3], 0, 3)
  n.of.mt.gene = sum(mt.gene.index)
  expr.binary.mt.gene.only = expr.binary[mt.gene.index,]
  unique.gene.counts.mt = colSums(expr.binary.mt.gene.only)
  
  gene.counts = colSums(expr)
  expr.mt.gene.only = expr[mt.gene.index,]
  gene.counts.mt = colSums(expr.mt.gene.only)
  
  
  
  filter.index.3 = gene.counts.mt <= gene.counts*0.1
  filter.index.3 = unique.gene.counts.mt <= unique.gene.counts*0.1
  
  filter.index.1 * filter.index.2 * filter.index.3
}
