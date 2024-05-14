setwd("/bioNSW/lixin/hyp/single_cell_ATAC/BCC")
rm(list=ls())
gc()
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(Seurat)
library(stringr)
library(Cairo)
library(ggplot2)
library(chromVAR)
library(chromVARmotifs)
library(ArchR)
#Set the thread
addArchRThreads(threads=22)
library(motifmatchr)
library(BiocParallel)
register(SerialParam())
set.seed(99)
genome <- BSgenome.Hsapiens.UCSC.hg19
projBCC1 <- loadArchRProject("Save-elncRNA-BCC")
elncRNAMatrixBCC <- getMatrixFromProject(
  ArchRProj = projBCC1,
  useMatrix = "elncRNAMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE
)
saveRDS(elncRNAMatrixBCC,"/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/peakmtx/elncRNAMatrixBCC.rds")
##TF activity
se <- readRDS("/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/peakmtx/elncRNAMatrixBCC.rds")
names(se@assays@data@listData) <- "counts"
se <- se[rowSums(se@assays@data@listData[["counts"]])>0,]
se <- addGCBias(se, genome = genome)#Computes GC content
data("human_pwms_v1")#in R package chromVARmotifs
#match Motifs
matches=matchMotifs(human_pwms_v1,rowRanges(se),genome = "BSgenome.Hsapiens.UCSC.hg19")
##Computes deviations in chromatin accessibility across sets of annotations
dev <- computeDeviations(object = se, annotations = matches)
saveRDS(dev, "E:\\单细胞\\pbs\\deviation.rds")
#compute variability
metadata(dev)$Variability <- computeVariability(dev)
plotVariability(metadata(dev)$Variability)
plot(sort(metadata(dev)$Variability[["variability"]],decreasing = T),
     xlab="Sorted TFs",ylab="Variability",main="1764 TFs Variability",pch=20)
abline(h = sort(metadata(dev)$Variability[["variability"]],decreasing = T)[250],v = 250, col = "gray36")
axis(side=1,at=250,labels=250, font.axis = 2)
#add se
metadata(dev)$SummarizedExperiment <- se
#add matches
metadata(dev)$motifMatches <- matches
saveRDS(dev, "E:\\单细胞\\pbs\\chromVAR-Summarized-Experiment.rds")
##250 most variable TFs across all scATAC-seq clusters
varTF_i <- sort(head(order(dev@metadata[["Variability"]][["variability"]], decreasing = TRUE), 250))
varTF_score=deviationScores(dev)[varTF_i,]
cluster_TF=NULL
for (i in 1:20) {
  cluster=rowMedians(varTF_score[,which(cell_barcode$Clusters==paste0("Cluster",i))],na.rm = T)
  cluster_TF<-cbind(cluster_TF,cluster)
}
colnames(cluster_TF) = cluster_name
rownames(cluster_TF) = dev@metadata[["Variability"]][["name"]][varTF_i]
library(pheatmap)
labels_row=row.names(cluster_TF)
labels_row[-which(row.names(cluster_TF) %in% c("TCF3","TCF21","TCF23","BATF","IRF4","NR3C1",
                                               "RUNX3","EOMES","TP53","FOS","SMARCC1","TBX21",
                                               "JUN","CEBPB","YY1","GATA1","NFKB1","STAT1","SOX8",
                                               "RARA","GATA1","NR4A1","TCF4","IRF6","GATA5","ERF"))]=""
pheatmap(cluster_TF,clustering_method = "mcquitty", show_rownames = T,
         color = colorRampPalette(c("#436eee", "white","#FF7F00", "#EE0000"))(30),
         cellwidth = 15, cellheight = 2,fontsize_col = 8,fontsize_row = 5,
         labels_row=labels_row, treeheight_row = 0, treeheight_col = 0,main = "250 most variable TFs deviation score")