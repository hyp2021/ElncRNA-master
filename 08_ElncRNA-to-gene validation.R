#################################################################
#validate elncRNA-gene interactions
###############################################################
#ABC enhancer-gene link
abc_dt <- fread("AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
keep_cols <- c(
  "chr", "start", "end", "name", "class", 
  "activity_base", "TargetGene", 
  #"TargetGeneTSS", "TargetGeneExpression", "TargetGenePromoterActivityQuantile", "TargetGeneIsExpressed", 
  "distance", 
  #"isSelfPromoter", "hic_contact", "powerlaw_contact", "powerlaw_contact_reference", "hic_contact_pl_scaled", 
  #"hic_pseudocount", "hic_contact_pl_scaled_adj", "ABC.Score.Numerator", "powerlaw.Score.Numerator", "powerlaw.Score", 
  "ABC.Score", 
  "CellType")
abc_dt <- abc_dt[,..keep_cols]
library(GenomicRanges)
library(rtracklayer)
abc_gr <- makeGRangesFromDataFrame(abc_dt, keep.extra.columns=TRUE, ignore.strand=TRUE, 
                                   seqnames.field="chr", start.field="start", end.field="end")
#enhancer lncRNA
lncRNA <- read.table(file = "noblacklist_enhancer_lncRNA.bed")
#lncRNA TSS
lncRNA_TSS <- GenomicRanges::GRanges(
  seqnames = lncRNA$V1,
  ranges = IRanges(lncRNA$V2, end = lncRNA$V3),
  strand = lncRNA$V4,
  symbol = lncRNA$V5
)
#lncRNA TSS overlap enhancer
o <- findOverlaps(lncRNA_TSS, abc_gr, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any")
lncRNA <- lncRNA[unique(queryHits(o)),]
#elnc_enhancer_gene
E2P2G <- abc_gr[subjectHits(o)]
E2P2G$elncrna <- lncRNA_TSS$symbol[queryHits(o)]
saveRDS(E2P2G,"ABC.rds")
scREG <- fread("scEnhancer_elncRNA_gene.txt")
names(scREG)[1] <- "elncrna"
names(scREG)[4] <- "celltype"
ABC <- data.frame(elncrna = E2P2G$elncrna,
                  enhancer = paste0(seqnames(E2P2G),":",ranges(E2P2G)),
                  gene = E2P2G$TargetGene,
                  ABC.Score = E2P2G$ABC.Score,
                  celltype = E2P2G$CellType)
ABC <- ABC[(E2P2G$distance<=250000)&(E2P2G$TargetGene %in% scREG$gene),]
ABC <- data.frame(ABC[!duplicated(ABC[,c(1,3,5)]),c(1,3,5)])
write.table(ABC,"ABC_elnc_gene.txt",row.names=F,col.names=T,quote=F,sep="\t",fileEncoding = "UTF-8")
length(unique(paste0(ABC$elncrna,ABC$gene,ABC$celltype)))
overlap <- length(intersect(unique(paste0(scREG$elncrna,scREG$gene)),unique(paste0(ABC$elncrna,ABC$gene))))
#total 250kb elnc-gene links
load("protein_coding_gene.RData")
genes <- GenomicRanges::GRanges(
  seqnames = protein_coding_gene$seqnames,
  ranges = IRanges(protein_coding_gene$start, end = protein_coding_gene$end),
  strand = protein_coding_gene$strand,
  gene_id = protein_coding_gene$gene_id,
  symbol = protein_coding_gene$gene_name
)
genes <- genes[(genes$symbol %in% scREG$gene),]
genes_TSS <- resize(resize(genes, 1), 2 * 250000 + 1, "center")
#enhancer lncRNA
#lncRNA TSS
lncRNA_TSS <- resize(resize(lncRNA_TSS, 1), 2 * 250000 + 1, "center")

#dist-match
o <- findOverlaps(lncRNA_TSS, genes_TSS, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any")
dm <- data.frame(elncrna=lncRNA_TSS[queryHits(o)]$symbol,
                 gene=genes[subjectHits(o)]$symbol)
dm <- dm[(dm$elncrna %in% unique(c(scREG$elncrna,ABC$elncrna)))&(dm$gene %in% unique(c(scREG$gene,ABC$gene))),]#46717
#hypergeometric
x<-overlap
m<-length(unique(paste0(scREG$elncrna,scREG$gene)))#17287
n<-nrow(dm)#46717
k<-length(unique(paste0(ABC$elncrna,ABC$gene)))#4512
SS=x/min(m,k)#0.7316046
p.value<-phyper(x-1,m,n-m,k,lower.tail=FALSE, log.p = T)#-1358.2952
#cell type separately
#ABC category (The 131 human cell types are classified into seven major cell type categories, such as mononuclear phagocytes, B cells, T cells, other haematopoietic cells, fibroblasts, epithelial cells, and other cell types.)
category <- fread("category.txt",header = F)
category <- category[(category$V1 != ""),]
nameidx <- grep("\\(",category$V1)
listname <- category$V1[nameidx]
category <- category[which(category$V1 %in% c(ABC$celltype,listname)),]
listname <- stringr::str_split(listname, pattern = " \\(", n = 2, simplify = TRUE)[,1]
gaps <- diff(c(nameidx,nrow(category)+1))-1
categorylist <- list()
for (i in 1:length(unique(listname))) {
  categorylist[[listname[i]]] <- category$V1[(nameidx[i]+1):(nameidx[i]+gaps[i])]
}
#single-cell elnc-gene links
celltypedetails <- readRDS("details_of_cell_type.rds")
names(celltypedetails)[1] <- "celltype"
names(celltypedetails)[2] <- "description"
categorylist2 <- list()
categorylist2[[listname[1]]] <- celltypedetails$celltype[setdiff(grep("CD14|THP|mono|Mono|DC|endritic|Macro|macro|phag|Antigen presenting cells|yeloid",celltypedetails$description),grep("CD14 negative NK",celltypedetails$description))]
categorylist2[[listname[2]]] <- celltypedetails$celltype[setdiff(grep("B |GM12878|lasma|Pre B|Memory$|Naive$|-B$",celltypedetails$description),grep("dendritic",celltypedetails$description))]
categorylist2[[listname[3]]] <- celltypedetails$celltype[setdiff(grep("T |Treg|Tfh|Th2|Thymocytes|CD4|CD8|Jurkat|K562",celltypedetails$description),grep("GIST",celltypedetails$description))]
categorylist2[[listname[4]]] <- celltypedetails$celltype[setdiff(grep("HSC|asophil|BM |blast cell|Progenitor|CLP.1|CMP|GMP|LMPP|MCP|MDP|CD34|MPAL|rythro|NK|kill|like|matopoi|progenitor|Neut|Platelet|pDC|spleen|egakaryocyte|Late.Eryth",celltypedetails$description),grep("CD14_negative_NK",celltypedetails$description))]
categorylist2[[listname[5]]] <- celltypedetails$celltype[setdiff(grep("fibro|Fibro|Astrocyte|HT1080",celltypedetails$description),grep("CD14_negative_NK",celltypedetails$description))]
categorylist2[[listname[6]]] <- celltypedetails$celltype[setdiff(grep("A549|alveolar|Basal cells|Ciliated |Club cells|Ureteric bud cells|eratinocyte|MCF7|MCF10A|T24|HeLa|pitheli|Tumor|cancer|GIST|Adrenocortical|intestine|pancreas|stomach|muscle",celltypedetails$description),grep("CD14_negative_NK",celltypedetails$description))]
categorylist2[[listname[7]]] <- celltypedetails$celltype[-which(celltypedetails$celltype %in% c(categorylist2[[listname[1]]],categorylist2[[listname[2]]],categorylist2[[listname[3]]],categorylist2[[listname[4]]],categorylist2[[listname[5]]],categorylist2[[listname[6]]]))]
for (i in 1:length(unique(listname))) {
  ABCi <- ABC[(ABC$celltype %in% categorylist[[i]]),]
  scREGi <- scREG[(scREG$celltype %in% categorylist2[[i]]),]
  a <- length(unique(paste0(ABCi$elncrna,ABCi$gene)))
  b <- length(unique(paste0(scREGi$elncrna,scREGi$gene)))
  overlap <- length(intersect(unique(paste0(scREGi$elncrna,scREGi$gene)),unique(paste0(ABCi$elncrna,ABCi$gene))))
  ss <- overlap/min(a,b)
  p.value<-phyper(overlap-1,b,46717-b,a,lower.tail=FALSE, log.p = T)
  print(paste0(listname[i]," : scREG : ",b," ABC : ",a," overlap : ",overlap," SS : ",ss," log(P) : ",p.value))
}