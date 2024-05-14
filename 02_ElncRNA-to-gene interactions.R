setwd("/bioNSW/lixin/hyp/single_cell_ATAC/scATAC")

##Run Cicero co-accessibility
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
library(Seurat)
library(dplyr)
rm(list=ls())
gc()
set.seed(102)

flank <- 250*10^3
corCutOff <- 0.1
#input each cell type
indata <- readRDS("/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/scATAC_mtx/celltype.rds")
#remove uncommon chromosome
chr <- paste0(c(paste0("chr",c(1:22,"X","Y"),"_","[0-9]"),paste0("chr",c(1:22,"X","Y"),":","[0-9]")), collapse="|")
indata <- indata[grep(chr, rownames(indata)),]

# binarize the matrix
indata@x[indata@x > 0] <- 1
#less than 200 peaks were removed
indata <- indata[,Matrix::colSums(indata)>=200]
#Ensure there are no peaks included with zero reads
indata <- indata[Matrix::rowSums(indata) != 0,] 
# format cell info
cellinfo <- data.frame(
  row.names = paste0(1:ncol(indata)),
  cells = colnames(indata)
)

# format peak info
peaks_bed <- stringr::str_split(paste0(rownames(indata)), pattern = "_|:|-", n = 3, simplify = TRUE)
rownames(peaks_bed) <- paste0(peaks_bed[,1],"_",peaks_bed[,2],"_",peaks_bed[,3])
peakinfo <- data.frame(
  row.names = rownames(peaks_bed),
  site_name = rownames(peaks_bed),
  chr = peaks_bed[,1],
  bp1 = peaks_bed[,2],
  bp2 = peaks_bed[,3]
)

rownames(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

#LSI
#Reduced Dimensions
#TF-IDF
message("Computing Term Frequency IDF...")
freqs <- t(t(indata)/Matrix::colSums(indata))
idf   <- as(log(1 + ncol(indata) / Matrix::rowSums(indata)), "sparseVector")
tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
#Calculate SVD then LSI
message("Computing SVD using irlba...")
svd <- irlba::irlba(tfidf, 50, 50)
svdDiag <- matrix(0, nrow=50, ncol=50)
diag(svdDiag) <- svd$d
matSVD <- t(svdDiag %*% t(svd$v))
rownames(matSVD) <- colnames(indata)
colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
#UMAP
umap_coords <- RunUMAP(object = matSVD, dims = 2:50)
umap_coords <- Embeddings(umap_coords)
# make CDS
cds <-  suppressWarnings(newCellDataSet(tfidf,
                                        phenoData = methods::new("AnnotatedDataFrame", data = cellinfo),
                                        featureData = methods::new("AnnotatedDataFrame", data = peakinfo),
                                        expressionFamily=negbinomial.size(),
                                        lowerDetectionLimit=0))

fData(cds)$bp1 <- as.numeric(fData(cds)$bp1)
fData(cds)$bp2 <- as.numeric(fData(cds)$bp2)
cds <- cds[order(fData(cds)$chr, fData(cds)$bp1),]
#detectGenes：Counts how many cells each feature in a CellDataSet object that are detectably expressed above a minimum threshold.
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
#cicero requires input of a CDS object
ciceroObj <- make_cicero_cds(cds, k = 50, reduced_coordinates = umap_coords)#Function to generate an aggregated input CDS for cicero, where k represents the data of how many cells need to be merged per bin

#Compute Cicero co-accessibility
data("human.hg19.genome")
sample_genome <- human.hg19.genome[1:24,]
interactions <- run_cicero(ciceroObj, sample_genome)
interactions=subset(interactions,coaccess>corCutOff)

saveRDS(interactions,paste0("/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/interaction/","celltype.rds"))

interactions <- rbind(interactions,data.frame(Peak1=interactions$Peak2,
                                              Peak2=interactions$Peak1,
                                              coaccess=interactions$coaccess))
i_enhancer.bed <- stringr::str_split(paste0(interactions[,1]), pattern = "_|:|-|[|]", n = 3, simplify = TRUE)
i_enhancer.bed <- data.frame(
  chr = i_enhancer.bed[,1],
  bp1 = i_enhancer.bed[,2],
  bp2 = i_enhancer.bed[,3],
  enhancer = paste0(i_enhancer.bed[,1],"_",i_enhancer.bed[,2],"_",i_enhancer.bed[,3])
)

i_promoter.bed <- stringr::str_split(paste0(interactions[,2]), pattern = "_|:|-|[|]", n = 3, simplify = TRUE)
i_promoter.bed <- data.frame(
  chr = i_promoter.bed[,1],
  bp1 = i_promoter.bed[,2],
  bp2 = i_promoter.bed[,3],
  promoter = paste0(i_promoter.bed[,1],"_",i_promoter.bed[,2],"_",i_promoter.bed[,3])
)
E_center = (as.numeric(i_enhancer.bed$bp1)+as.numeric(i_enhancer.bed$bp2))/2
P_center = (as.numeric(i_promoter.bed$bp1)+as.numeric(i_promoter.bed$bp2))/2
idx <- which(abs(E_center-P_center)<=250*10^3)
interactions <- data.frame(
  Peak1 = i_enhancer.bed$enhancer,
  Peak2 = i_promoter.bed$promoter,
  coaccess = interactions[,3]
)
interactions <- interactions[idx,]
interactions <- distinct(interactions)
i_enhancer.bed <- distinct(i_enhancer.bed)
i_promoter.bed <- distinct(i_promoter.bed)
saveRDS(interactions,paste0("/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/interaction_E_P/","celltype.rds"))
write.table(i_enhancer.bed,file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/scEnhancer_allpromoter.bed",row.names=F,col.names=F,quote=F,sep="\t")
write.table(i_promoter.bed,file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/scEnhancer_allenhancer.bed",row.names=F,col.names=F,quote=F,sep="\t")

#bedtools intersect -a noblacklist_enhancer_lncRNA.bed -b scEnhancer_allenhancer.bed -wa -wb > scEnhancer_elncRNA_peak.bed
#bedtools intersect -a promoter.bed -b scEnhancer_allpromoter.bed -wa -wb > scEnhancer_promoter_peak.bed

setwd("/bioNSW/lixin/hyp/single_cell_ATAC/scATAC")
library(dplyr)
#annotate gene
gene_peak=read.table(file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/scEnhancer_promoter_peak.bed")
gene_peak=data.frame(site_name=gene_peak$V8,
                     gene=gene_peak$V4)
gene_peak=data.frame(distinct(gene_peak))
write.table(unique(gene_peak$gene),file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/gene.txt",row.names=F,col.names=F,quote=F,sep="\t")
#annotate enhancer lncRNA
noblacklist_enhancer_lncRNA_peaks=read.table(file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/scEnhancer_elncRNA_peak.bed")
noblacklist_enhancer_lncRNA_peaks = data.frame(site_name=noblacklist_enhancer_lncRNA_peaks$V11,
                                               gene=noblacklist_enhancer_lncRNA_peaks$V5)
noblacklist_enhancer_lncRNA_peaks <- data.frame(distinct(noblacklist_enhancer_lncRNA_peaks))
write.table(unique(noblacklist_enhancer_lncRNA_peaks$gene),
            file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/elncRNA.txt",row.names=F,col.names=F,quote=F,sep="\t")
#immune gene
immune_gene=read.table("/bioNSW/lixin/hyp/single_cell_ATAC/bed/GeneList.txt", header = T, sep = "\t", quote = "", comment.char = "#")
write.table(unique(immune_gene$Symbol),file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/immune_gene.txt",row.names=F,col.names=F,quote=F,sep="\t")
#annotate peak
annotate_peak=rbind(gene_peak,noblacklist_enhancer_lncRNA_peaks)

#elncRNA-gene co-accessibility
connections <- interactions
connections[] <- lapply(connections, function(x) {
  if (is.factor(x)) x <- as.character(x)
  return(x)
})

interactions=left_join(connections, noblacklist_enhancer_lncRNA_peaks, by =c("Peak1" = "site_name"))
interactions=left_join(interactions,gene_peak,by =c("Peak2" = "site_name"))
interactions=subset(interactions,gene.x!="NA"&gene.y!="NA",c(4,5,3))
interactions=data.frame(elncRNA=interactions$gene.x,
                        gene=interactions$gene.y,
                        coaccess=interactions$coaccess)
interactions[] <- lapply(interactions, function(x) {
  if (is.factor(x)) x <- as.character(x)
  return(x)
})
interactions <- interactions[interactions[,1]!=interactions[,2],]
interactions=interactions %>% group_by(elncRNA,gene) %>% mutate(coaccess = median(coaccess)) %>% ungroup()
interactions <- data.frame(interactions[!duplicated(interactions[,c(1:3)]),])
#filter immune_gene_elncRNA interactions
immune_interactions=subset(interactions,gene %in% immune_gene$Symbol)
write.table(immune_interactions,file = paste0("elncRNA_immune_interactions/","celltype.txt"),row.names=F,col.names=T,quote=F,sep="\t")
#filter all_gene_elncRNA interactions
write.table(interactions,file = paste0("elncRNA_allgene_interactions/","celltype.txt"),row.names=F,col.names=T,quote=F,sep="\t")