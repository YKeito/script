"TF_TargetGenes_Dendrogram.R"
#install.packages("factoextra")
#package----
library(dplyr)
library(tidyr)
library(factoextra)
library(ggplot2)
library(stringr)
#Input data----
AtRegNet.Table <- readRDS(file = "~/nishioka/RDS/AtRegNetTable.rds")
MasterTable <- readRDS(file = "~/nishioka/RDS/MasterTable.rds")
BRM.GSEV <- readRDS(file = "~/nishioka/RDS/BRM_GSEV.rds")
#Focus on TF-DEG target genes----
T.Venn <- readRDS(file = "~/nishioka/RDS/Nishioka_VennDiagram.rds")
T.AGI <- T.Venn$GeneList
AtRegNet.Table <- AtRegNet.Table %>% filter(AGI %in% T.AGI)
#processing data----
df <- AtRegNet.Table %>% select(annotation, MCLNum, pvalue)
BRM.GSEV <- BRM.GSEV %>% select(MCLNum, pvalue)
BRM.GSEV <- data.frame(annotation = rep("BRM", times = nrow(BRM.GSEV))) %>% cbind(BRM.GSEV)
df <- df %>% rbind(BRM.GSEV)
#significantly enriched target genes in cluster----
T.data <- df %>% 
  filter(pvalue < 0.05)
T.MCLNum <- T.data %>% 
  distinct(MCLNum) %>% 
  unlist(use.names = F)
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 70]
df <- df %>% 
  filter(MCLNum %in% T.MCLNum)
BRM.MCLNum <- df %>% 
  filter(annotation == "BRM", pvalue < 0.05) %>%
  distinct(MCLNum) %>% unlist(use.names = F)
BRM.GRN <- df %>% 
  filter(MCLNum %in% BRM.MCLNum, pvalue < 0.05) %>% 
  left_join(y = AtRegNet.Table %>% select("AGI", "annotation", "MCLNum", "AGRIS", "N", "M", "n", "x"), 
            by = c("annotation" = "annotation", "MCLNum" = "MCLNum")
            )
#Hierarchical Clustering----
T.matrix <- df %>% spread(key = MCLNum, value = pvalue)
rownames(T.matrix) <- T.matrix$annotation
T.matrix <- T.matrix %>% select(-annotation)
T.logic <- T.matrix < 0.05
T.logic <- T.logic %>% apply(MARGIN = 1, FUN = sum)
T.matrix <- T.matrix[T.logic != 0, ]
T.annotation <- rownames(T.matrix)
df <- df %>% 
  filter(annotation %in% T.annotation)
#BRM data
df <- df %>% 
  filter(MCLNum %in% BRM.MCLNum)
res <- hclust(dist(T.matrix), method = "ward.D2")
df$annotation <- factor(x = df$annotation,
                        levels = res$labels[res$order],
                        ordered = TRUE)
res.MCL <- hclust(dist(t(T.matrix)), method = "ward.D2")
df$MCLNum <- factor(x = df$MCLNum,
                    levels = res.MCL$labels[res.MCL$order],
                    ordered = TRUE)
#Divide p-value into five----
df$pvalue[df$pvalue > 5e-2] <- 1
df$pvalue[df$pvalue <= 5e-2 & df$pvalue > 5e-5] <- 0.67
df$pvalue[df$pvalue <= 5e-5 & df$pvalue > 5e-10] <- 0.34
df$pvalue[df$pvalue <= 5e-10] <- 0
#GSEV:plotting heatmap----
g <- ggplot(data = df, aes(x = MCLNum, y = annotation, fill = pvalue))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient(high="white",low="red", na.value = "white")
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g <- g + theme(axis.ticks.y = element_blank())
g <- g + theme(legend.position = "top")
g <- g + theme(legend.position = 'none')
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + theme(legend.title = element_blank())
g <- g + theme(axis.text = element_text(size=20))
g <- g + coord_flip()
ggsave(filename = "~/nishioka/Image/GSEV/BRM&AtRegNet_TFDEGs_GSEV.png", plot = g, width = 20, height = 8)
write.table(x = T.matrix, file = "~/nishioka/Image/BRM&AGRIS_GSEVScore.txt", sep = "\t", quote = F, row.names = F)
write.table(x = df, file = "~/nishioka/Table/BRM&AGRIS_GSEVTable.txt", sep = "\t", quote = F, row.names = F)
#dendrogram----
#Hclust; annotation
res.hc <- eclust(x = T.matrix,
                 "hclust",
                 k = ncol(T.matrix),
                 method = "euclidean",
                 graph = FALSE
)
#Dendrogram----
#annotation
g <- fviz_dend(res.hc,
               cex = 0.5,
               color_labels_by_k = TRUE,
               show_labels = TRUE,
               ggtheme = theme_classic(),
               horiz = FALSE,
               #rect = TRUE,
               rect_fill = TRUE,
               type = "rectangle",
               main = NULL
)
g <- g + theme(axis.text=element_text(size=20))
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
T.HCL <- data.frame(annotation = names(res.hc$cluster),
                    Hcluster = res.hc$cluster,
                    order = res.hc$order,
                    row.names = NULL)
#option----
#k:クラスター数の指定
#k_colors, palette:色指定
#show_labels:横軸の名前を表示するか
#color_labels_by_k:クラスターのグループに応じて色を付けるか
#horiz縦軸と横軸を反転するかしないか
#rect:クラスターを囲む
#rect_fill:rectで囲んだ範囲を色塗する
#type:図の形式の指定
#lwd:dendrogramの線の太さ変更
#label_cols:下の文字の色指定
#save----
ggsave(filename = "~/nishioka/Image/Dendrogram/BRM&AtRegNet_DEGs_Dendrogram.png", g, width = 20, height = 10)
saveRDS(T.HCL, "~/nishioka/RDS/BRM&AtRegNet_DEGs_Dendrogram.rds")
write.table(T.HCL, "~/nishioka/Table/BRM&AtRegNet_DEGs_Dendrogram.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = T)
write.table(BRM.GRN, "~/nishioka/Table/BRM&AtRegNet_DEGs_GRN.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = T)