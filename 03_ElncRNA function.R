frefro <- c("HSC_fresh_BM0828","HSC_frozen_BM0828")
frefro_fdir <- paste0(frefro, "_interaction.txt")
frefro_inter <- NULL
for (i in 1:length(frefro_fdir)) {
  inter <- read.table(paste0("E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\elncRNA_allgene_interactions\\", frefro_fdir[i]),header = T)
  inter <- data.frame(inter,cell = frefro[i])
  frefro_inter <- rbind(frefro_inter,inter)
}
write.table(frefro_inter,"E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\results\\frefro\\frefro_inter.txt",row.names=F,col.names=T,quote=F,sep="\t")
#interact enrichment
frefro_inter <- read.table("E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\results\\frefro\\frefro_inter.txt",header = T,sep="\t")
library(org.Hs.eg.db)
library(clusterProfiler)

#Enrichment of frefro intersect
interactions <- frefro_inter[duplicated(frefro_inter[,1:2]),]
#elncRNA co-accessibility gene
lncRNA2gene=unique(interactions$gene)

geneid <-  AnnotationDbi::select(org.Hs.eg.db, keys=lncRNA2gene, keytype="SYMBOL",
                                 columns="ENTREZID")
entrezID=as.integer(geneid$ENTREZID)

#KEGG enrichment
KEGG <- enrichKEGG(gene = entrezID,
                   organism = "hsa",
                   keyType = "kegg",
                   pAdjustMethod = "none",
                   pvalueCutoff = 1,
                   qvalueCutoff  = 1)
KEGG <- setReadable(KEGG,  org.Hs.eg.db, keyType="ENTREZID")
KEGG_result=KEGG[KEGG$pvalue<0.01]
write.table(KEGG_result,paste0("E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\results\\frefro\\","frefro"," EnrichmentKEGG.txt"), row.names =FALSE, quote = F, sep = "\t")
#Draw a bubble chart 
p=barplot(
  KEGG,
  x = "GeneRatio",
  color = "pvalue",
  showCategory = nrow(KEGG_result),
  size = NULL,
  split = NULL,
  font.size = 9,
  title = paste0("HSC"," EnrichmentKEGG"),
  label_format = 30,
  orderBy = "x"
)
print(p)
KEGG <- enrichGO(gene = entrezID,
                 org.Hs.eg.db,
                 ont = "BP",
                 pAdjustMethod = "none",
                 pvalueCutoff = 1,
                 qvalueCutoff  = 1)
KEGG <- setReadable(KEGG,  org.Hs.eg.db, keyType="ENTREZID")
KEGG_result=KEGG[KEGG$pvalue<0.01]
write.table(KEGG_result,paste0("E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\results\\frefro\\","frefro"," EnrichmentGO.txt"), row.names =FALSE, quote = F, sep = "\t")
#Draw a bubble chart 
p=barplot(
  KEGG,
  x = "GeneRatio",
  color = "pvalue",
  showCategory = 5,
  size = NULL,
  split = NULL,
  font.size = 9,
  title = paste0("HSC"," EnrichmentGO"),
  label_format = 30,
  orderBy = "x"
)
print(p)
#single elncRNA
#MIR210HG
elnc <- c("MIR210HG","CTD-2126E3.3","HCG25","RP11-1398P2.1","ERVK3-1","LINC00963")
res <- NULL
for (i in 1:length(unique(elnc))) {
  lncRNA2gene=unique(interactions$gene[which(interactions$elncRNA %in% elnc[i])])
  geneid <-  AnnotationDbi::select(org.Hs.eg.db, keys=lncRNA2gene, keytype="SYMBOL",
                                   columns="ENTREZID")
  entrezID=as.integer(geneid$ENTREZID)
  KEGG <- enrichGO(gene = entrezID,
                   org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "none",
                   pvalueCutoff = 1,
                   qvalueCutoff  = 1)
  KEGG <- setReadable(KEGG,  org.Hs.eg.db, keyType="ENTREZID")
  KEGG_result=KEGG[KEGG$pvalue<0.01]
  KEGG_result$Cluster <- elnc[i]
  res <- rbind(res,KEGG_result)
  write.table(res,paste0("E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\results\\frefro\\","elnc"," EnrichmentGO.txt"), row.names =FALSE, quote = F, sep = "\t")
}
library(dplyr)
enrich <- res %>% 
  group_by(Cluster) %>% 
  top_n(n = 3, wt = -pvalue) %>% 
  filter(Cluster %in% elnc)

dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = dt$Description)
library(ggplot2)
#theme
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(),
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)

p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =c('#6bb9d2', '#d55640','#ff8787','#7ed9e0','#ff9999','#00ccff')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "-Log10(pvalue)", y = paste0(rev(names(table(dt$Cluster))),collapse="    "), title = "ElncRNA-associated function") +
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) +
  geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust = 2.5, color=rep(c('#6bb9d2', '#d55640','#ff8787','#7ed9e0','#ff9999','#00ccff'),times=table(dt$Cluster))) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme

p