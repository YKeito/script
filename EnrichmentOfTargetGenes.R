"~/bigdata/yasue/nishioka/script/EnrichmentOfTargetGenes.R"
#load package----
library(dplyr)
library(tidyr)
library(ggplot2)
#Input data----
AtRegNet.Table <- readRDS(file = "~/bigdata/yasue/nishioka/RDS/AtRegNetTable.rds")
MasterTable <- readRDS(file = "~/bigdata/yasue/nishioka/RDS/MasterTable.rds")
DofList <- read.table(file = "~/bigdata/yasue/nishioka/base/allDOFTF.txt", sep = "\t", header = T, stringsAsFactors = F)
HDZIPIII <- read.table(file = "~/bigdata/yasue/nishioka/base/HDZIPIII.txt", sep = "\t", header = T, stringsAsFactors = F)
#processing data----
T.logic <- MasterTable %>% select(X0_2d_FDR005:X0_10d_FDR005) == "Yes"
T.logic <- T.logic %>% apply(MARGIN = 1, FUN = sum) != 0
T.DEGs <- MasterTable$AGI[T.logic]
Dof.AGI <- DofList$AGI %>% intersect(y = T.DEGs)
HDZIPIII.AGI <- HDZIPIII$AGI %>% intersect(y = T.DEGs)
T.AGI <- c(Dof.AGI, HDZIPIII.AGI)
AtDofRegNet.Table <- AtRegNet.Table %>% filter(AGI %in% T.AGI)
AtDofRegNet.Table <- AtDofRegNet.Table %>% filter(M >= 10, n >= 10)
df <- AtDofRegNet.Table %>% select(annotation, MCLNum, pvalue)
T.matrix <- df %>% spread(key = MCLNum, value = pvalue)
rownames(T.matrix) <- T.matrix$annotation
T.matrix <- T.matrix %>% select(-annotation)
res <- hclust(dist(T.matrix), method = "ward.D2")
df$annotation <- factor(x = df$annotation,
                        levels = res$labels[res$order],
                        ordered = TRUE)
res.MCL <- hclust(dist(t(T.matrix)), method = "ward.D2")
df$MCLNum <- factor(x = df$MCLNum,
                    levels = res.MCL$labels[res.MCL$order],
                    ordered = TRUE)
df$pvalue[df$pvalue > 5e-2] <- 1
df$pvalue[df$pvalue <= 5e-2 & df$pvalue > 5e-5] <- 0.67
df$pvalue[df$pvalue <= 5e-5 & df$pvalue > 5e-10] <- 0.34
df$pvalue[df$pvalue <= 5e-10] <- 0
T.MCLNum <- df$MCLNum %>% unique()
T.data <- c()
i <- 1
for(i in i:length(T.MCLNum)){
  test <- sum(df$pvalue[df$MCLNum == T.MCLNum[i]] == 1) != df$annotation %>% unique() %>% length()
  if(test){
    T.data <- rbind(T.data,
                    df %>% filter(MCLNum == T.MCLNum[i])
                    )
  }
  i <- i+1
}
g <- ggplot(data = T.data, aes(x = MCLNum, y = annotation, fill = pvalue))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient(high="white",low="red", na.value = "white")
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g <- g + theme(axis.ticks.y = element_blank())
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + theme(legend.title = element_blank())
g <- g + theme(legend.text = element_text(size=7))
plot(g)
ggsave(filename = "~/bigdata/yasue/nishioka/Image/DOFHDZIPIII_EnrichmentOfTargetGenes.png", plot = g)
write.table(x = T.matrix, file = "~/bigdata/yasue/nishioka/Table/DEGsDofHDZIPIII_GSEVScore.txt", sep = "\t", quote = F, row.names = F)
write.table(x = AtDofRegNet.Table, 
            file = "~/bigdata/yasue/nishioka/Table/DEGsDofHDZIPIII_GSEVTable.txt", sep = "\t", quote = F, row.names = F)