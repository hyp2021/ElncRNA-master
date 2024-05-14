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
######################
#ALL
#####################
#Set the required gene and genome annotations
blacklist <- import.bed("/bioNSW/lixin/hyp/single_cell_ATAC/bed/consensusBlacklist.bed")
genomeAnnotation <- createGenomeAnnotation(
  genome = BSgenome.Hsapiens.UCSC.hg19,
  chromSizes = NULL,
  blacklist = blacklist,
  filter = TRUE,
  filterChr = c("chrM")
)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                       OrgDb = org.Hs.eg.db)

#Record fragments' path in a vector and pass it as a parameter to createArrowFiles()
fdir <- list.files("BCCFragments", pattern = "_fragments.tsv.gz$")
inputFiles <- paste0("BCCFragments/", fdir)
names(inputFiles) <- str_extract(fdir,"(?<=_).+?(?=_fragments.tsv.gz)")
##
# reformatFragmentFiles(
#   fragmentFiles = c("GSM3722035_CD4_Naive_fragments.tsv.gz","GSM3722034_CD4_Memory_fragments.tsv.gz"),
#   checkChrPrefix = getArchRChrPrefix()
# )
#Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 4,
  minFrags = 1000,
  maxFrags = 1e+20,
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = F
)

#Creating An ArchRProject
projBCC1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Save-scATAC-ALL-Granja-2019",
  copyArrows = TRUE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)
#add blacklist
projBCC1@genomeAnnotation$blacklist <- blacklist
#add ENCODE GeneScoreMatrix
#import gencode.v19.annotation.gtf_withproteinids
hg19_annotation= rtracklayer::import("/bioNSW/lixin/hyp/single_cell_ATAC/enhancerlnc/gencode.v19.annotation.gtf.gz",format = "gtf")

hg19_annotation=as.data.frame(hg19_annotation)
#select protein coding gene
protein_coding_gene=subset(hg19_annotation,(type=="gene")&(gene_type=="protein_coding"),c("seqnames","start","end","strand","gene_id","gene_name"))
save(protein_coding_gene,file = "/bioNSW/lixin/hyp/single_cell_ATAC/enhancerlnc/protein_coding_gene.RData")

genes <- GenomicRanges::GRanges(
  seqnames = protein_coding_gene$seqnames,
  ranges = IRanges(protein_coding_gene$start, end = protein_coding_gene$end),
  strand = protein_coding_gene$strand,
  gene_id = protein_coding_gene$gene_id,
  symbol = protein_coding_gene$gene_name
)
blacklist <- import.bed("/bioNSW/lixin/hyp/single_cell_ATAC/bed/consensusBlacklist.bed")
projBCC1 <- addGeneScoreMatrix(
  input = projBCC1,
  genes = genes,
  matrixName = "GeneScoreMatrix",
  excludeChr = c("chrY", "chrM"),
  blacklist = blacklist,
  force = T,
  logFile = createLogFile("addGeneScoreMatrix")
)
saveArchRProject(ArchRProj = projBCC1, outputDirectory = "Save-scATAC-ALL-Granja-2019", load = FALSE)

projBCC1 <- loadArchRProject("Save-scATAC-ALL-Granja-2019")

#subset 37818 known cell
matrix_dir = "/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/peakmtx/"
barcode.path <- paste0(matrix_dir, "GSE129785_scATAC-TME-All.cell_barcodes.txt")
features.path <- paste0(matrix_dir, "GSE129785_scATAC-TME-All.peaks.txt")
cell_barcode = read.delim(barcode.path,
                          header = TRUE,
                          stringsAsFactors = FALSE)

projBCC1 <- loadArchRProject("Save-scATAC-ALL-Granja-2019")
projBCC1 <- subsetArchRProject(
  ArchRProj = projBCC1,
  cells = cell_barcode$Group_Barcode,
  outputDirectory = "Save-elncRNA-BCC",
  dropCells = TRUE,
  force = TRUE
)

projBCC1 <- loadArchRProject("Save-elncRNA-BCC")

projBCC1 <- addCellColData(ArchRProj = projBCC1, data = cell_barcode$Clusters,
                           cells = cell_barcode$Group_Barcode, name = "Clusters")
Cluster_names = c("1-Naive CD4 T","2-Th17","3-Tfh","4-Treg","5-Naive CD8 T",
                  "6-Th1","7-Memory CD8 T","8-CD8 TEx","9-Effector CD8 T","10-NK1",
                  "11-NK2","12-B","13-Plasma B","14-Myeloid cells","15-Endothelial",
                  "16-Fibroblasts","17-Tumor 1","18-Tumor 2","19-Tumor 3","20-Tumor 4")
projBCC1$Cluster_names <- mapLabels(projBCC1$Clusters, newLabels = Cluster_names, oldLabels = paste0("Cluster",1:20))


#add 580789 peak
feature.names = read.delim(features.path,
                           header = TRUE,
                           stringsAsFactors = FALSE)
#GRanges
chr=t(as.data.frame(strsplit(feature.names$Feature,"_")))
peaks <- GenomicRanges::GRanges(seqnames = as.vector(chr[,1]),
                                ranges = IRanges::IRanges(start = as.numeric(chr[,2]),
                                                          width = 501))
projBCC1 <- addPeakSet(
  ArchRProj = projBCC1,
  peakSet = peaks,
  force = T
)
#add peak matrix
projBCC1 <- addPeakMatrix(projBCC1,
                          ceiling = 10^9,
                          binarize = FALSE)


#samples reside in which clusters
cM <- confusionMatrix(paste0(projBCC1$Cluster_names), paste0(projBCC1$Sample))
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p

set.seed(1)
#Iterative Latent Semantic Indexing
projBCC1 <- addIterativeLSI(
  ArchRProj = projBCC1,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  iterations = 4,
  scaleDims = T,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
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
  totalFeatures = 580789,
  filterQuantile = 0.995,
  binarize = TRUE,
  force = T
)

#UMAP 2:50
projBCC1 <- addUMAP(
  ArchRProj = projBCC1, 
  reducedDims = "IterativeLSI", 
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
#Identifying Background Peaks
projBCC1 <- addBgdPeaks(projBCC1,force = TRUE)
projBCC1 <- addDeviationsMatrix(
  ArchRProj = projBCC1, 
  peakAnnotation = "Motif",
  force = TRUE
)

#Computing Impute Weights Using Magic -> marker UMAP
projBCC1 <- addImputeWeights(projBCC1)

#Making Pseudo-bulk Replicates
projBCC1<- addGroupCoverages(ArchRProj = projBCC1, groupBy = "Cluster_names",force = T)

#trajectory

#elncRNA

#add elncRNA matrix
#bedtools intersect -a noblacklist_enhancer_lncRNA.bed -b peaks_bed_BCC.bed -wa -wb > noblacklist_enhancer_lncRNA_peaks_BCC.bed
noblacklist_enhancer_lncRNA_peaks=read.table("/bioNSW/lixin/hyp/single_cell_ATAC/bed/noblacklist_enhancer_lncRNA_peaks_BCC.bed")
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
  force = TRUE
)

saveArchRProject(ArchRProj = projBCC1, outputDirectory = "Save-elncRNA-BCC", load = FALSE)

projBCC1 <- loadArchRProject("Save-elncRNA-BCC")
markersPeaks <- getMarkerFeatures(
  ArchRProj = projBCC1, 
  useMatrix = "elncRNAMatrix", 
  groupBy = "Cluster_names",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 10^9
)
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  nLabel = 2,
  nPrint = 2,
  transpose = F
)
plotPDF(heatmapPeaks, name = "elncRNA-Marker-Heatmap", width = 8, height = 6, ArchRProj = projBCC1, addDOC = FALSE)
saveRDS(getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1"),"Save-elncRNA-BCC/Plots/elncRNAmarkerBCC.rds")