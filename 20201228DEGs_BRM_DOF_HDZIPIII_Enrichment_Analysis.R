#load package----
library(tidyverse)
#input data----
DEGs.Ave <- read.table("~/IMO/base/20201217_Aveppm_DEGs.txt", sep = "\t", header = T, quote = "")
AtRegNet.Table <- readRDS(file = "~/IMO/RDS/20201222AtRegNetTable.rds")
T.DEGs <- DEGs.Ave$AGI
DEGs.AtRegNet.Table <- AtRegNet.Table %>% filter(AGI %in% T.DEGs)
DEGs.AtRegNet.Table <- DEGs.AtRegNet.Table %>% select(- AtRegNet_TFFamily)
BRM.GSEV <- readRDS(file = "~/IMO/RDS/20201222BRM_GSEV.rds")
BRM.GSEV <- BRM.GSEV %>% mutate(AGI = "AT2G46020", annotation = "BRM")
T.data <- DEGs.AtRegNet.Table %>% rbind(BRM.GSEV)
Check.list <- read.table(file = "~/IMO/base/DOF_iMO.txt", sep = "\t", header = T, quote = "")
#processing data----
T.AGI <- c(Check.list$AGI, "AT2G46020")
T.data <- T.data %>% 
  filter(M >= 10, n >= 10) %>% 
  select(- Intersection_AGI) %>% 
  mutate(count = ifelse(pvalue <= 0.05, 1, 0)) %>% 
  filter(AGI %in% T.AGI)

#T.MCLNum <- T.data %>% 
#  group_by(MCLNum) %>% 
#  summarise(Summary = sum(count)) %>% 
#  filter(Summary == 0) %>% 
#  select(MCLNum) %>% 
#  unlist(use.names = F)

#T.AGI <- T.data %>% 
#  group_by(AGI) %>% 
#  summarise(Summary = sum(count)) %>% 
#  filter(Summary == 0) %>% 
#  select(AGI) %>% 
#  unlist(use.names = F)


T.AGI.summary <- T.data %>% 
  group_by(AGI) %>% 
  summarise(Summary = sum(count)) %>% 
  filter(Summary != 0)

T.MCL.summary <- T.data %>% 
  group_by(MCLNum) %>% 
  summarise(Summary = sum(count)) %>% 
  filter(Summary != 0)


T.data <- T.data %>% 
  filter(AGI %in% T.summary$AGI, MCLNum %in% T.MCL.summary$MCLNum)

df <- T.data %>% select(annotation, MCLNum, pvalue)
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
g <- g + coord_flip()
plot(g)
ggsave(filename = "~/IMO/Image/20201228DEGs_BRM_Dof_HSZIPIII_GSEV.png", plot = g, width = 1.5, height = 6)
write.table(x = T.matrix, file = "~/IMO/Table/20201228DEGs_BRM_Dof_HSZIPIII_GSEVScore.txt", sep = "\t", quote = F, row.names = F)
write.table(x = DEGs.AtRegNet.Table, 
            file = "~/IMO/Table/20201228DEGs_BRM_Dof_HSZIPIII_GSEVTable.txt", sep = "\t", quote = F, row.names = F)