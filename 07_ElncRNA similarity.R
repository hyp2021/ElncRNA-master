elncRNA_cell <- read.table("E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\results\\elncRNA_cell.txt", header = T, sep = "\t", quote = "", comment.char = "#")
elncRNA_cell=data.frame(elncRNA=elncRNA_cell[,1],
                        cell=elncRNA_cell[,4])

Jmat <- matrix(nrow = length(unique(elncRNA_cell$cell)),ncol = length(unique(elncRNA_cell$cell)))
rownames(Jmat) <- unique(elncRNA_cell$cell)
colnames(Jmat) <- unique(elncRNA_cell$cell)
for (i in 1:length(unique(elncRNA_cell$cell))) {
  for (j in 1:length(unique(elncRNA_cell$cell))) {
    Jmat[i,j] <- length(intersect(subset(elncRNA_cell,cell==rownames(Jmat)[i])$elncRNA,
                                  subset(elncRNA_cell,cell==colnames(Jmat)[j])$elncRNA
    ))/length(union(subset(elncRNA_cell,cell==rownames(Jmat)[i])$elncRNA,
                    subset(elncRNA_cell,cell==colnames(Jmat)[j])$elncRNA))
  }
}
saveRDS(Jmat,"cell_Jaccard_similarity.rds")
write.table(Jmat,file = "ElncRNA similarity.txt",row.names=T,col.names=T,quote=F,sep="\t")
#532 cell types
library(pheatmap)
set.seed(1)
annotation_col <- data.frame(
  Tissue = information$Tissue,
  Batch = information$batch,
  Disease =information$Health_state
)
rownames(annotation_col) <- information$Cell.type
pdf("E:\\houyaopan\\scATAC-seq\\scElncRNA_result\\results\\aapaper\\cor532.pdf",width=20, height=15)
pheatmap(
  Jmat,
  clustering_method = "ward.D",
  annotation_col = annotation_col,border_color = NA,treeheight_row = 0, treeheight_col = 0
)
dev.off()