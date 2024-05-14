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
set.seed(1)
addArchRThreads(threads=22)
##################################
#Hematopoiesis
###############################
#subset 61806 cell
matrix_dir = "/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/peakmtx/"
barcode.path <- paste0(matrix_dir, "GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt")
features.path <- paste0(matrix_dir, "GSE129785_scATAC-Hematopoiesis-All.peaks.txt")
cell_barcode = read.delim(barcode.path,
                          header = TRUE,
                          stringsAsFactors = FALSE)

idxPass <- which(cell_barcode$Group_Barcode %in% projBCC1$cellNames)

projBCC1 <- subsetArchRProject(
  ArchRProj = projBCC1,
  cells = cell_barcode$Group_Barcode[idxPass],
  outputDirectory = "Save-elncRNA-Hematopoiesis",
  dropCells = TRUE,
  force = TRUE
)
projBCC1 <- loadArchRProject("Save-elncRNA-Hematopoiesis")
#cell cluster 2019
idxPass <- which(cell_barcode$Group_Barcode %in% projBCC1$cellNames)
projBCC1 <- addCellColData(ArchRProj = projBCC1, data = cell_barcode$Clusters[idxPass],
                           cells = cell_barcode$Group_Barcode[idxPass], name = "Clusters")
Cluster_names = c("1-HSC/MPP",
                  "2-MEP",
                  "3-CMP/BMP",
                  "4-LMPP",
                  "5-CLP",
                  "6-Pro-B",
                  "7-Pre-B",
                  "8-GMP",
                  "9-MDP",
                  "10-pDC",
                  "11-cDC",
                  "12-Monocyte 1",
                  "13-Monocyte 2",
                  "14-Naive B",
                  "15-Memory B",
                  "16-Plasma cell",
                  "17-Basophil",
                  "18-Immature NK",
                  "19-Mature NK1",
                  "20-Mature NK2",
                  "21-Naive CD4 T1",
                  "22-Naive CD4 T2",
                  "23-Naive Treg",
                  "24-Memory CD4 T",
                  "25-Treg",
                  "26-Naive CD8 T1",
                  "27-Naive CD8 T2",
                  "28-Naive CD8 T3",
                  "29-Central memory CD8 T",
                  "30-Effector memory CD8 T",
                  "31-Gamma delta T")
projBCC1$Cluster_names <- mapLabels(projBCC1$Clusters, newLabels = Cluster_names, oldLabels = paste0("Cluster",1:31))

#peak
feature.names = read.delim(features.path,
                           header = TRUE,
                           stringsAsFactors = FALSE)
#GRanges
chr <- t(as.data.frame(strsplit(feature.names$Feature,"_")))
peaks <- GenomicRanges::GRanges(seqnames = as.vector(chr[,1]),
                                ranges = IRanges::IRanges(start = as.numeric(chr[,2]),
                                                          width = 501))
projBCC1 <- addPeakSet(
  ArchRProj = projBCC1,
  peakSet = peaks,
  force = FALSE
)
#add peak matrix
projBCC1 <- addPeakMatrix(projBCC1,
                          ceiling = 10^9,
                          binarize = FALSE)

set.seed(1)
#Iterative Latent Semantic Indexing
projBCC1 <- addIterativeLSI(
  ArchRProj = projBCC1,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  iterations = 4,
  scaleDims = T,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.05), 
    sampleCells = 10000,
    maxClusters = NULL,
    n.start = 10
  ),
  UMAPParams = list(n_neighbors =30,
                    min_dist = 0.3,
                    metric = "cosine",
                    verbose =FALSE),
  varFeatures = 25000, 
  dimsToUse = 1:50,
  LSIMethod = 1,
  corCutOff = 0.75,
  totalFeatures = 571400,
  filterQuantile = 0.995,
  binarize = TRUE,
  force = T
)
projBCC1 <- addHarmony(
  ArchRProj = projBCC1,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:50,
  scaleDims = T,
  corCutOff = 0.75,
  name = "Harmony",
  groupBy = "Sample",
  force = T,
  do_pca = F,
  lambda = 5,
  sigma = 0.01,
  max.iter.harmony = 1
)
#UMAP
projBCC1 <- addUMAP(
  ArchRProj = projBCC1, 
  reducedDims = "Harmony", 
  name = "UMAP", 
  dimsToUse = 2:50,
  nNeighbors = 30, 
  minDist = 0.3, 
  metric = "cosine",
  force = T
)
p1 <- plotEmbedding(ArchRProj = projBCC1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projBCC1, colorBy = "cellColData", name = "Cluster_names", embedding = "UMAP")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projBCC1, addDOC = FALSE, width = 5, height = 5)

#peak motif analysis
#Motif Annotations
projBCC1 <- addMotifAnnotations(ArchRProj = projBCC1, motifSet = "cisbp", name = "Motif",force = TRUE)

#chromVAR
if("Motif" %ni% names(projBCC1@peakAnnotation)){
  projBCC1 <- addMotifAnnotations(ArchRProj = projBCC1, motifSet = "cisbp", name = "Motif")
}
projBCC1 <- addBgdPeaks(projBCC1,force = TRUE)
projBCC1 <- addDeviationsMatrix(
  ArchRProj = projBCC1, 
  peakAnnotation = "Motif",
  force = TRUE
)

#addImputeWeight
projBCC1 <- addImputeWeights(projBCC1)

#addGroupCoverages
projBCC1<- addGroupCoverages(ArchRProj = projBCC1, groupBy = "Cluster_names",force = TRUE)

#trajectory

#elncRNA

#add elncRNA matrix
noblacklist_enhancer_lncRNA_peaks=read.table("/bioNSW/lixin/hyp/single_cell_ATAC/bed/noblacklist_enhancer_lncRNA_peaks_Hematopoiesis.bed")
elncRNApeaks <- GenomicRanges::GRanges(seqnames = noblacklist_enhancer_lncRNA_peaks$V1,
                                       ranges = IRanges::IRanges(start = as.numeric(noblacklist_enhancer_lncRNA_peaks$V10),
                                                                 width = 501))

mcols(elncRNApeaks)$name <- paste0(noblacklist_enhancer_lncRNA_peaks$V8)
projBCC1 <- addFeatureMatrix(
  input = projBCC1,
  features = elncRNApeaks,
  matrixName = "elncRNAMatrix",
  ceiling = 10^9,
  binarize = FALSE,
  force = TRUE,
  logFile = createLogFile("addFeatureMatrix")
)

saveArchRProject(ArchRProj = projBCC1, outputDirectory = "Save-elncRNA-Hematopoiesis", load = FALSE)
#########################################################################
#Integrative pseudo-time analyses for elncRNA-gene
#########################################################################
##B traj
set.seed(1)
addArchRThreads(threads=22)
projBCC1 <- loadArchRProject("Save-elncRNA-Hematopoiesis")
maxDist = 250000
corCutOff = 0.35
varCutOff = 0.75
fdrCutoff = 0.001

trajectory <- c("1-HSC/MPP","4-LMPP","5-CLP","6-Pro-B","7-Pre-B","14-Naive B","15-Memory B","16-Plasma cell")
projBCC1 <- addTrajectory(
  ArchRProj = projBCC1, 
  name = "B_traj", 
  groupBy = "Cluster_names",
  trajectory = trajectory, 
  embedding = "UMAP",
  force = TRUE
)
p <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "cellColData", name = "B_traj")
plotPDF(p, name = "Plot-B_traj-Traj-UMAP.pdf", ArchRProj = projBCC1, addDOC = FALSE, width = 5, height = 5)
trajEM  <- getTrajectory(ArchRProj = projBCC1, name = "B_traj", useMatrix = "elncRNAMatrix", log2Norm = TRUE)
trajGSM  <- getTrajectory(ArchRProj = projBCC1, name = "B_traj", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
corEM_GSM <- correlateTrajectories(
  seTrajectory1 = trajEM,
  seTrajectory2 = trajGSM,
  corCutOff = corCutOff,
  varCutOff1 = varCutOff,
  varCutOff2 = varCutOff,
  removeFromName1 = NULL,
  removeFromName2 = NULL,
  useRanges = TRUE,
  fix1 = "center",
  fix2 = "start",
  maxDist = maxDist,
  log2Norm1 = TRUE,
  log2Norm2 = TRUE,
  force = FALSE,
  logFile = createLogFile("correlateTrajectories")
)
#subset our corresponding trajectory
corEM_GSM[[1]] <- subset(corEM_GSM[[1]], FDR < fdrCutoff)

df <- data.frame(corEM_GSM[[1]],
                 name1 = rownames(trajEM)[corEM_GSM[[1]]$idx1],
                 name2 = rownames(trajGSM)[corEM_GSM[[1]]$idx2])

saveRDS(df,"/bioNSW/lixin/hyp/single_cell_ATAC/BCC/Save-elncRNA-Hematopoiesis/Plots/corEM_GSM_B_traj.rds")
#label top 10
labelMarkers <- df[unique(c(grep("CXCR4|HOXA9|MEF2C|EBF1|PAX5|RUNX1|IL7R|RAG2|PRDM1|SPIB|IKZF|STAT5|CD19|LYN|CD34|CD20|BCL11A|SPI1|IRF8|TCF3|AIOLOS",df$name2))),]

trajEM2 <- trajEM[df$idx1, ]
trajGSM2 <- trajGSM[df$idx2, ]

trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajEM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajEM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, 
                             rowOrder = rowOrder,labelMarkers = rownames(trajEM)[labelMarkers$idx1],labelTop = 0)
ht2 <- plotTrajectoryHeatmap(trajGSM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0,
                             rowOrder = rowOrder,labelMarkers = rownames(trajGSM)[labelMarkers$idx2],labelTop = 0)
plotPDF(ht1 + ht2, name = "Plot-B_traj-Traj-Heatmaps.pdf", ArchRProj = projBCC1, addDOC = FALSE, width = 10, height = 8)
saveRDS(ht2,"/bioNSW/lixin/hyp/single_cell_ATAC/BCC/Save-elncRNA-Hematopoiesis/Plots/trajGSMht2.rds")
p1 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "GeneScoreMatrix", name = c("RUNX1"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "horizonExtra", log2Norm = TRUE)
p2 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "GeneScoreMatrix", name = c("HOXA9"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "horizonExtra", log2Norm = TRUE)
p3 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "GeneScoreMatrix", name = c("PAX5"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "horizonExtra", log2Norm = TRUE)
p4 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "GeneScoreMatrix", name = c("BCL11A"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "horizonExtra", log2Norm = TRUE)
p5 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "GeneScoreMatrix", name = c("LYN"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "horizonExtra")
p <- ggAlignPlots(p1[[2]],p2[[2]],p3[[2]],p4[[2]],p5[[2]], type = "v", draw = F)
plotPDF(p, name = "Plot-B_traj-Traj-gene.pdf", ArchRProj = projBCC1, addDOC = FALSE, width = 8, height = 20)

p1 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "elncRNAMatrix", name = c("4637_AP000688.29"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "blueYellow", log2Norm = TRUE)
p2 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "elncRNAMatrix", name = c("1927_HOXA11-AS"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "blueYellow", log2Norm = TRUE)
p3 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "elncRNAMatrix", name = c("2372_RP11-297B17.3"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "blueYellow", log2Norm = TRUE)
p4 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "elncRNAMatrix", name = c("5492_AC007381.3"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "blueYellow", log2Norm = TRUE)
p5 <- plotTrajectory(projBCC1, trajectory = "B_traj", colorBy = "elncRNAMatrix", name = c("2190_RP11-318K15.2"),imputeWeights = getImputeWeights(projBCC1), continuousSet = "blueYellow", log2Norm = TRUE)
p <- ggAlignPlots(p1[[2]],p2[[2]],p3[[2]],p4[[2]],p5[[2]], type = "v", draw = F)
plotPDF(p, name = "Plot-B_traj-Traj-elncRNA.pdf", ArchRProj = projBCC1, addDOC = FALSE, width = 8, height = 20)

saveArchRProject(ArchRProj = projBCC1, outputDirectory = "Save-elncRNA-Hematopoiesis", load = FALSE)