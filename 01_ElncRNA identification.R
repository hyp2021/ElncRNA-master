#elncRNA
library(rtracklayer)
#import lncRNA
lncRNA=rtracklayer::import("gencode.v19.long_noncoding_RNAs.gtf")
lncRNA=as.data.frame(lncRNA)
#select transcript
lncRNA_transcript=subset(lncRNA,type=="transcript",c("seqnames","start","end","strand","transcript_id","gene_name"))
save(lncRNA_transcript,file = "lncRNA_transcript.RData")
#import gencode.v19.annotation.gtf_withproteinids
hg19_annotation= rtracklayer::import("gencode.v19.annotation.gtf_withproteinids",format = "gtf")
hg19_annotation=as.data.frame(hg19_annotation)
#select protein-coding gene
protein_coding_gene=subset(hg19_annotation,(type=="gene")&(gene_type=="protein_coding"),c("seqnames","start","end","strand","gene_id","gene_name"))
save(protein_coding_gene,file = "protein_coding_gene.RData")
#promoter(tssWindow=2kb)
TSS1=subset(protein_coding_gene,strand=="+",c("seqnames","start","gene_name"))
colnames(TSS1)=c("seqnames","TSS","gene_name")
TSS2=subset(protein_coding_gene,strand=="-",c("seqnames","end","gene_name"))
colnames(TSS2)=c("seqnames","TSS","gene_name")
TSS=rbind(TSS1,TSS2)
TSS[] <- lapply(TSS, function(x) {
  if (is.factor(x)) x <- as.character(x)
  return(x)
})
promoter=data.frame(seqnames=TSS$seqnames,
                    start=TSS$TSS-2000,
                    end=TSS$TSS+2000,
                    gene_name=TSS$gene_name)
write.table(promoter,file = "bed\\promoter.bed",row.names=F,col.names=F,quote=F,sep="\t")
setwd("/bioNSW/lixin/hyp/single_cell_ATAC/BCC")
promoter=read.table(file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/promoter.bed")
write.table(unique(promoter$V4),
            file = "/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/scElncRNA_result/gene.txt",row.names=F,col.names=F,quote=F,sep="\t")
#lncRNA bed
lncRNA_bed=lncRNA_transcript[,c("seqnames","start","end","strand","gene_name")]
write.table(lncRNA_bed,file = "bed\\lncRNA.bed",row.names=F,col.names=F,quote=F,sep="\t")
#non-promoter lncRNA
#bedtools intersect -a lncRNA.bed -b promoter.bed -v > non_promoter_lncRNA.bed
#enhancer
path <- "K27ac_bed\\"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x)}) 
temp<-data.frame()
for (i in 1:length(fileNames)){
  temp<-rbind(temp,subset(data[[i]],V4=="7_EnhG1"|V4=="8_EnhG2"|V4=="9_EnhA1"|V4=="10_EnhA2"))
}
write.table(temp,file = "enhancer.bed",row.names=F,col.names=F,quote=F,sep="\t")
#enhancer lncRNA
non_promoter_lncRNA=read.table("bed\\non_promoter_lncRNA.bed")
library(dplyr)
non_promoter_lncRNA=distinct(non_promoter_lncRNA)
lncRNA_TSS1=subset(non_promoter_lncRNA,V4=="+",c(V1,V2,V2,V4,V5,V2,V3))
colnames(lncRNA_TSS1)=c("seqnames","TSS1","TSS2","strand","gene_name","start","end")
lncRNA_TSS2=subset(non_promoter_lncRNA,V4=="-",c(V1,V3,V3,V4,V5,V2,V3))
colnames(lncRNA_TSS2)=c("seqnames","TSS1","TSS2","strand","gene_name","start","end")
#lncRNA TSS
lncRNA_TSS=rbind(lncRNA_TSS1,lncRNA_TSS2)
write.table(lncRNA_TSS,file = "bed\\lncRNA_TSS.bed",row.names=F,col.names=F,quote=F,sep="\t")
#lncRNA TSS overlap enhancer
#bedtools intersect -a lncRNA_TSS.bed -b enhancer.bed -wa > enhancer_lncRNA.bed
enhancer_lncRNA=read.table("bed\\enhancer_lncRNA.bed")
enhancer_lncRNA=distinct(enhancer_lncRNA)
enhancer_lncRNA=data.frame(seqnames=enhancer_lncRNA$V1,
                           start=enhancer_lncRNA$V6,
                           end=enhancer_lncRNA$V7,
                           strand=enhancer_lncRNA$V4,
                           gene_name=enhancer_lncRNA$V5,
                           TSS1=enhancer_lncRNA$V2,
                           TSS2=enhancer_lncRNA$V3)
write.table(enhancer_lncRNA,file = "bed\\enhancer_lncRNA.bed",row.names=F,col.names=F,quote=F,sep="\t")
#enhancer lncRNA remove blacklist 
#bedtools intersect -a enhancer_lncRNA.bed -b consensusBlacklist.bed -v > noblacklist_enhancer_lncRNA.bed
noblacklist_enhancer_lncRNA=read.table(file = "/bioNSW/lixin/hyp/single_cell_ATAC/bed/noblacklist_enhancer_lncRNA.bed")
write.table(unique(noblacklist_enhancer_lncRNA$V5),
            file = "/bioNSW/lixin/hyp/single_cell_ATAC/scATAC/scElncRNA_result/elncRNA.txt",row.names=F,col.names=F,quote=F,sep="\t")
noblacklist_enhancer_lncRNA=read.table("bed\\noblacklist_enhancer_lncRNA.bed")
noblacklist_enhancer_lncRNA=distinct(noblacklist_enhancer_lncRNA)
noblacklist_enhancer_lncRNA=data.frame(seqnames=noblacklist_enhancer_lncRNA$V1,
                                       TSS1=noblacklist_enhancer_lncRNA$V6,
                                       TSS2=noblacklist_enhancer_lncRNA$V7,
                                       strand=noblacklist_enhancer_lncRNA$V4,
                                       gene_name=noblacklist_enhancer_lncRNA$V5,
                                       start=noblacklist_enhancer_lncRNA$V2,
                                       end=noblacklist_enhancer_lncRNA$V3)
write.table(noblacklist_enhancer_lncRNA,file = "bed\\noblacklist_enhancer_lncRNA.bed",row.names=F,col.names=F,quote=F,sep="\t")
noblacklist_enhancer_lncRNA=read.table(file = "noblacklist_enhancer_lncRNA.bed")
noblacklist_enhancer_lncRNA=cbind(noblacklist_enhancer_lncRNA,paste0(1:9652,"_",noblacklist_enhancer_lncRNA[,5]))
write.table(noblacklist_enhancer_lncRNA,file = "noblacklist_enhancer_lncRNA.bed",row.names=F,col.names=F,quote=F,sep="\t")