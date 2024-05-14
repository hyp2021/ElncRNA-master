#################################################################
#validate elncRNA-gene interactions
###############################################################
#Enhancer atlas 2.0 enhancer-gene link
fdir <- list.files("atlas2_EP", pattern = "_EP.txt$")
inputFiles <- paste0("atlas2_EP/", fdir)
names(fdir) <- gsub("_EP.txt","",fdir)
atlas2EP <- NULL
for(i in 1:length(fdir)){
  interactions=read.table(paste0("atlas2_EP/", fdir[i]),header = F)
  subgene=stringr::str_split(paste0(interactions[,1]), pattern = "\\$", n = 5, simplify = TRUE)
  subenhancer=stringr::str_split(paste0(subgene[,1]), pattern = ":|-|_", n = 4, simplify = TRUE)
  bed=data.frame(seqnames=subenhancer[,1],start=subenhancer[,2],end=subenhancer[,3],gene=subgene[,2],score=interactions[,2],cell=names(fdir)[i],distance=abs((as.numeric(subenhancer[,2])+as.numeric(subenhancer[,3]))/2-as.numeric(subgene[,4])))
  bed=bed[which(bed$distance<=250*10^3),]
  atlas2EP <- rbind(atlas2EP,bed)
}
write.table(atlas2EP,file = paste0("atlas2EP",".bed"),row.names=F,col.names=F,quote=F,sep="\t")
#bedtools intersect -a noblacklist_enhancer_lncRNA.bed -b atlas2EP.bed -wa -wb > atlas2_elncRNA_gene.bed
atlas2EG=read.table(paste0("atlas2_elncRNA_gene.bed"),header = F)
atlas2EG=data.frame(elncRNA=atlas2EG$V5,
                    gene=atlas2EG$V12,
                    score=atlas2EG$V13,
                    celltissue=atlas2EG$V14)
gene_peak=read.table(file = "bed/promoter.bed")
atlas2EG=atlas2EG[which(atlas2EG$gene %in% gene_peak$V4),]
library(dplyr)
atlas2EG <- atlas2EG %>% group_by(celltissue,elncRNA,gene) %>% mutate(score = median(score)) %>% ungroup()
atlas2EG <- data.frame(atlas2EG[!duplicated(atlas2EG[,c(1,2,4,3)]),])
write.table(atlas2EG,"atlas2_elncRNA_gene.txt",row.names=F,col.names=T,quote=F,sep="\t")
#scEnhancer elncRNA-gene links
immune <- information$Cell.type
immune_fdir <- paste0(immune, "_interaction.txt")
immune_inter <- NULL
for (i in 1:length(immune_fdir)) {
  inter <- read.table(paste0("elncRNA_allgene_interactions/", immune_fdir[i]),header = TRUE)
  inter <- data.frame(inter,
                      cell = immune[i],
                      tissue = information$Tissue[match(immune[i],information$Cell.type)])
  immune_inter <- rbind(immune_inter,inter)
}
write.table(immune_inter,"scEnhancer_elncRNA_gene.txt",row.names=F,col.names=T,quote=F,sep="\t")

library(pROC)
#Heart
atlas2EG=read.table(paste0("atlas2_elncRNA_gene.txt"),header = TRUE)
a2=atlas2EG[grep("heart",atlas2EG$celltissue),]
sc=immune_inter[grep("Heart",immune_inter$tissue),]

sc=data.frame(links=paste0(sc$elncRNA,sc$gene),
              score=sc$coaccess)
sc <- sc %>% 
  dplyr::group_by(links) %>% 
  dplyr::summarise(score = median(score)
  ) %>% ungroup()
a2=data.frame(links=paste0(a2$elncRNA,a2$gene),
              score=a2$score)
a2 <- a2 %>% 
  dplyr::group_by(links) %>% 
  dplyr::summarise(score = median(score)
  ) %>% ungroup()
#ROC
dataroc=data.frame(links=sc$links)
dataroc$a2=0
dataroc$sc=0
dataroc$a2[which(dataroc$links %in% a2$links)]=1
dataroc$sc[match(sc$links,dataroc$links)]=sc$score

roc1 <- roc(dataroc$a2,dataroc$sc, levels = c(0, 1), direction = "<")
plot(roc1,col='red',legacy.axes=T)
#Jaccard
intersects=data.frame(Heart=c(length(sc$links),length(a2$links),length(intersect(sc$links,a2$links))))
#Kidney
a2=atlas2EG[grep("kidney",atlas2EG$celltissue),]
sc=immune_inter[grep("Kidney",immune_inter$tissue),]
roc2 <- roc(dataroc$a2,dataroc$sc, levels = c(0, 1), direction = "<")
plot(roc2,add=T,col='blue')
intersects$Kidney=c(length(sc$links),length(a2$links),length(intersect(sc$links,a2$links)))
#placenta
a2=atlas2EG[grep("placenta",atlas2EG$celltissue),]
sc=immune_inter[grep("Placenta",immune_inter$tissue),]
roc3 <- roc(dataroc$a2,dataroc$sc, levels = c(0, 1), direction = "<")
plot(roc3,add=T,col='brown')
intersects$Placenta=c(length(sc$links),length(a2$links),length(intersect(sc$links,a2$links)))
#Spleen
a2=atlas2EG[grep("Spleen",atlas2EG$celltissue),]
sc=immune_inter[grep("Spleen",immune_inter$tissue),]
roc4 <- roc(dataroc$a2,dataroc$sc, levels = c(0, 1), direction = "<")
plot(roc4,add=T,col='green')
intersects$Spleen=c(length(sc$links),length(a2$links),length(intersect(sc$links,a2$links)))
#spinal
a2=atlas2EG[grep("spinal|Osteoblast",atlas2EG$celltissue),]
sc=immune_inter[grep("Bone Marrow",immune_inter$tissue),]
roc5 <- roc(dataroc$a2,dataroc$sc, levels = c(0, 1), direction = "<")
plot(roc5,add=T,col='pink')
intersects$Bone=c(length(sc$links),length(a2$links),length(intersect(sc$links,a2$links)))
legend("bottomright",legend = c(auc(roc1),auc(roc2),auc(roc3),auc(roc4),auc(roc5)),col = c("red","blue","brown","green","pink"),lty = 1)

#NC Hi-C enhancer-gene links
GM12878=read.table(paste0("ceRNA/lib25_GM12878.enh_1e-2.bed"),header = F)
GM12878$cell="GM12878"
K562=read.table(paste0("ceRNA/lib25_K562.enh_1e-2.bed"),header = F)
K562$cell="K562"
interactions=rbind(GM12878,K562)
subgene=stringr::str_split(paste0(interactions[,4]), pattern = ":|Kb", n = 3, simplify = TRUE)
bed=data.frame(seqnames=interactions[,1],start=interactions[,2],end=interactions[,3],gene=subgene[,1],cell=interactions$cell,distance=abs(as.numeric(subgene[,2])))
bed=bed[which(bed$distance<=250),]
write.table(bed,file = paste0("GM12878_K562",".bed"),row.names=F,col.names=F,quote=F,sep="\t")
#bedtools intersect -a noblacklist_enhancer_lncRNA.bed -b GM12878_K562.bed -wa -wb > GM12878_K562_elncRNA_gene.bed
atlas2EG=read.table(paste0("GM12878_K562_elncRNA_gene.bed"),header = F)
atlas2EG=data.frame(elncRNA=atlas2EG$V5,
                    gene=atlas2EG$V12,
                    celltissue=atlas2EG$V13)
gene_peak=read.table(file = "bed/promoter.bed")
atlas2EG=atlas2EG[which(atlas2EG$gene %in% gene_peak$V4),]
library(dplyr)
atlas2EG <- data.frame(atlas2EG[!duplicated(atlas2EG[,c(1,2,3)]),])
write.table(atlas2EG,"GM12878_K562_elncRNA_gene.txt",row.names=F,col.names=T,quote=F,sep="\t")

#K562
a2=atlas2EG[grep("K562",atlas2EG$celltissue),1:2]
sc=immune_inter[grep("K562",immune_inter$cell),1:2]
sc <- data.frame(sc[!duplicated(sc[,c(1,2)]),])
scu=paste0(sc$elncRNA,sc$gene)
a2u=paste0(a2$elncRNA,a2$gene)
#Simpson
Simpson=length(intersect(scu,a2u))/min(length(scu),length(a2u))
#P10000
a<-0
for (n in 1:10000) {
  sc_random=paste0(sc$elncRNA,sample(sc$gene,length(sc$gene)))
  Simpson_random <- length(intersect(sc_random,a2u))/min(length(scu),length(a2u))
  if(Simpson_random > Simpson){a=a+1}
}
p=a/10000
Simpson10000p=data.frame(K562=c(length(scu),length(a2u),length(intersect(scu,a2u)),Simpson,p))
#GM12878
a2=atlas2EG[grep("GM12878",atlas2EG$celltissue),1:2]
sc=immune_inter[grep("GM12878",immune_inter$cell),1:2]
sc <- data.frame(sc[!duplicated(sc[,c(1,2)]),])
scu=paste0(sc$elncRNA,sc$gene)
a2u=paste0(a2$elncRNA,a2$gene)
#Simpson
Simpson=length(intersect(scu,a2u))/min(length(scu),length(a2u))
#P10000
a<-0
for (n in 1:10000) {
  sc_random=paste0(sc$elncRNA,sample(sc$gene,length(sc$gene)))
  Simpson_random <- length(intersect(sc_random,a2u))/min(length(scu),length(a2u))
  if(Simpson_random > Simpson){a=a+1}
}
p=a/10000
Simpson10000p$GM12878=c(length(scu),length(a2u),length(intersect(scu,a2u)),Simpson,p)