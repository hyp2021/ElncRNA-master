<<<<<<< Updated upstream
inter <- read.table(paste0("scEnhancer_elncRNA_gene.txt"), header = T, sep = "\t", quote = "", comment.char = "#")
information <- readRDS("details_of_cell_type.rds")
inter$health_state <- information$Health_state[match(inter$cell,information$Cell.type)]
cancer_selnc <- unique(setdiff(inter$elncRNA[-grep("Normal",inter$health_state)],inter$elncRNA[grep("Normal",inter$health_state)]))#133
normal <- unique(inter$elncRNA[grep("Normal",inter$health_state)])#2968
overlap <- unique(intersect(inter$elncRNA[-grep("Normal",inter$health_state)],inter$elncRNA[grep("Normal",inter$health_state)]))#1522
df=data.frame(table(paste0(inter$cell,"_",inter$elncRNA)))
df$elnc <- gsub(paste0(paste0(unique(information$Cell.type),"_"),collapse = "|"),"",df$Var1)
df$num.genes.Freq <- df$Freq
df$Freq[as.numeric(df$num.genes.Freq) <= 2] <- "1-2"
df$Freq[as.numeric(df$num.genes.Freq) > 2] <- "3+"
df$category[(df$elnc %in% cancer_selnc)] <- "Cancer-specific elncRNAs"
df$category[(df$elnc %in% normal)] <- "Normal elncRNAs"

cross_tab <- table(df$Freq, df$category)
cross_tab_df <- as.data.frame(cross_tab)
cross_tab_df$Proportion <- c(cross_tab_df$Freq[1]/sum(cross_tab_df$Freq[1:2])*100,cross_tab_df$Freq[2]/sum(cross_tab_df$Freq[1:2])*100,
                             cross_tab_df$Freq[3]/sum(cross_tab_df$Freq[3:4])*100,cross_tab_df$Freq[4]/sum(cross_tab_df$Freq[3:4])*100)

# Cancer specific elncRNAs link to less genes on average with statistical significance:
test <- wilcox.test(df$num.genes.Freq[(df$category %in% "Cancer-specific elncRNAs")],df$num.genes.Freq[(df$category %in% "Normal elncRNAs")],correct = F)
print(test)
print(test$p.value)
print(paste0("Cancer-specific elncRNAs link to less genes on average (",mean(df$num.genes.Freq[(df$category %in% "Cancer-specific elncRNAs")])," v. ",mean(df$num.genes.Freq[(df$category %in% "Normal elncRNAs")])," genes, p=",test$p.value,")"))

# Proportion of 3+ elncRNAs is less in cancer-specific elncRNAs relative to normal elncRNAs:
#col1: cancer, row1: 3+
res <- fisher.test(matrix(c(cross_tab_df$Freq[2], cross_tab_df$Freq[4], 
                            cross_tab_df$Freq[1], cross_tab_df$Freq[3]),
                          ncol = 2))
print(res)
print(res$p.value)
print(paste0("The proportion of cancer-specific elncRNAs linking to more genes is significantly less relative to the normal elncRNAs (p=",res$p.value,")"))
library(ggplot2)
ggplot(cross_tab_df,aes(x=Var1,y=Proportion,fill=Var2))+
  geom_bar(stat="identity",position = position_dodge(width=0.75),width=0.6)+
  theme_classic()+
  geom_text(aes(label=paste0(round(Proportion,1),"%")), vjust=-0.4, size=3.5,position = position_dodge(width=0.75))+
  scale_y_continuous(expand = c(0,0),limits=c(0,100))+scale_fill_manual(values = c("Cancer-specific elncRNAs" = "brown", "Normal elncRNAs" = "grey"))+
  labs(title = paste0("Proportion of elncRNAs per number of linked genes (p=",format(res$p.value,digits = 3),")"), x = "", fill = "")
#Mean number of linked genes per elncRNA
ggplot(df, aes(x=category, y=num.genes.Freq, fill = category)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",width=0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=0.15)+
  theme_classic()+  coord_cartesian(ylim=c(0.00,3.0))+
  scale_y_continuous(expand = c(0,0))+scale_fill_manual(values = c("Cancer-specific elncRNAs" = "brown", "Normal elncRNAs" = "grey"))+
  labs(title = paste0("Mean number of linked genes per elncRNA (p=",format(test$p.value, scientific = TRUE,digits = 3),")"), x = "", fill = "")
=======
inter <- read.table(paste0("scEnhancer_elncRNA_gene.txt"), header = T, sep = "\t", quote = "", comment.char = "#")
information <- readRDS("details_of_cell_type.rds")
inter$health_state <- information$Health_state[match(inter$cell,information$Cell.type)]
cancer_selnc <- unique(setdiff(inter$elncRNA[-grep("Normal",inter$health_state)],inter$elncRNA[grep("Normal",inter$health_state)]))#133
normal <- unique(inter$elncRNA[grep("Normal",inter$health_state)])#2968
overlap <- unique(intersect(inter$elncRNA[-grep("Normal",inter$health_state)],inter$elncRNA[grep("Normal",inter$health_state)]))#1522
df=data.frame(table(paste0(inter$cell,"_",inter$elncRNA)))
df$elnc <- gsub(paste0(paste0(unique(information$Cell.type),"_"),collapse = "|"),"",df$Var1)
df$num.genes.Freq <- df$Freq
df$Freq[as.numeric(df$num.genes.Freq) <= 2] <- "1-2"
df$Freq[as.numeric(df$num.genes.Freq) > 2] <- "3+"
df$category[(df$elnc %in% cancer_selnc)] <- "Cancer-specific elncRNAs"
df$category[(df$elnc %in% normal)] <- "Normal elncRNAs"

cross_tab <- table(df$Freq, df$category)
cross_tab_df <- as.data.frame(cross_tab)
cross_tab_df$Proportion <- c(cross_tab_df$Freq[1]/sum(cross_tab_df$Freq[1:2])*100,cross_tab_df$Freq[2]/sum(cross_tab_df$Freq[1:2])*100,
                             cross_tab_df$Freq[3]/sum(cross_tab_df$Freq[3:4])*100,cross_tab_df$Freq[4]/sum(cross_tab_df$Freq[3:4])*100)

# Cancer specific elncRNAs link to less genes on average with statistical significance:
test <- wilcox.test(df$num.genes.Freq[(df$category %in% "Cancer-specific elncRNAs")],df$num.genes.Freq[(df$category %in% "Normal elncRNAs")],correct = F)
print(test)
print(test$p.value)
print(paste0("Cancer-specific elncRNAs link to less genes on average (",mean(df$num.genes.Freq[(df$category %in% "Cancer-specific elncRNAs")])," v. ",mean(df$num.genes.Freq[(df$category %in% "Normal elncRNAs")])," genes, p=",test$p.value,")"))

# Proportion of 3+ elncRNAs is less in cancer-specific elncRNAs relative to normal elncRNAs:
#col1: cancer, row1: 3+
res <- fisher.test(matrix(c(cross_tab_df$Freq[2], cross_tab_df$Freq[4], 
                            cross_tab_df$Freq[1], cross_tab_df$Freq[3]),
                          ncol = 2))
print(res)
print(res$p.value)
print(paste0("The proportion of cancer-specific elncRNAs linking to more genes is significantly less relative to the normal elncRNAs (p=",res$p.value,")"))
library(ggplot2)
ggplot(cross_tab_df,aes(x=Var1,y=Proportion,fill=Var2))+
  geom_bar(stat="identity",position = position_dodge(width=0.75),width=0.6)+
  theme_classic()+
  geom_text(aes(label=paste0(round(Proportion,1),"%")), vjust=-0.4, size=3.5,position = position_dodge(width=0.75))+
  scale_y_continuous(expand = c(0,0),limits=c(0,100))+scale_fill_manual(values = c("Cancer-specific elncRNAs" = "brown", "Normal elncRNAs" = "grey"))+
  labs(title = paste0("Proportion of elncRNAs per number of linked genes (p=",format(res$p.value,digits = 3),")"), x = "", fill = "")
#Mean number of linked genes per elncRNA
ggplot(df, aes(x=category, y=num.genes.Freq, fill = category)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",width=0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=0.15)+
  theme_classic()+  coord_cartesian(ylim=c(0.00,3.0))+
  scale_y_continuous(expand = c(0,0))+scale_fill_manual(values = c("Cancer-specific elncRNAs" = "brown", "Normal elncRNAs" = "grey"))+
  labs(title = paste0("Mean number of linked genes per elncRNA (p=",format(test$p.value, scientific = TRUE,digits = 3),")"), x = "", fill = "")
>>>>>>> Stashed changes
