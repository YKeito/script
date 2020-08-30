"~/nishioka/script/DEGs_Enrichment_Analysis.R"
before <- proc.time()
#load package----
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
#Input data----
Nishioka.data <- read.table(file = "~/nishioka/base/20200317_log2FC&FDR.txt", sep = "\t", header = T, stringsAsFactors = F)
MasterTable <- readRDS(file = "~/nishioka/RDS/MasterTable.rds")
TAIR10 <- read.table(file = "~/nishioka/PublishedData/TAIR10.txt", stringsAsFactors = F, fill = T, sep = "\t", header = T)
TAIR10$AGI <- TAIR10$AGI %>% str_sub(end = -3)
TAIR10 <- TAIR10 %>% distinct(AGI)
#detecting DEGs----
DEGs <- list()
T.names <- c("2d", "4d", "7d", "10d")
df <- c()
i <- 1
for(i in i:length(T.names)){
  Nishioka.DEGs <- Nishioka.data %>% select(contains(T.names[i]))
  FDR.data <- Nishioka.DEGs %>% select(ends_with("qvalue"))
  FDR.logic <- FDR.data < 0.05
  FDR.logic <- FDR.logic %>% apply(MARGIN = 1, FUN = sum) == 2
  DEGs <- c(DEGs, list(Nishioka.data$reference_id[FDR.logic]))
  log2FC.data <- Nishioka.DEGs %>% select(ends_with("log2FC"))
  Up.logic <- log2FC.data > 0 & FDR.logic
  Down.logic <- log2FC.data < 0 & FDR.logic
  df <- rbind(df,
              data.frame(Sample = T.names[i],
                         Group = c("02Down", "01Up"),
                         value = c(sum(Down.logic), sum(Up.logic)))
  )
  i <- i+1
}
names(DEGs) <- T.names
#Enrichment analysis for TF Target genes----
T.MCLNum <- MasterTable$MCLNum[!is.na(MasterTable$MCLNum)] %>% unique()
N <- TAIR10 %>% nrow()
DEGs.Enrichment <- c()
i <- 1
for(i in i:length(DEGs)){
  Target.AGI <- DEGs[[i]]
  M <- Target.AGI %>% length()
  T.pvalue <- c()
  T.AGI <- c()
  n <- c()
  x <- c()
  j <- 1
  for(j in j:length(T.MCLNum)){
    T.data <- MasterTable %>% filter(MCLNum == T.MCLNum[j])
    n <- c(n, nrow(T.data))
    Cluster.AGI <- T.data$AGI
    T.AGI <- c(T.AGI, intersect(Target.AGI, Cluster.AGI))
    x <- c(x, length(intersect(Target.AGI, Cluster.AGI)))
    T.pvalue <- c(T.pvalue, phyper(x[j]-1, M, N-M, n[j], lower.tail = F))
    j <- j+1
  }
  DEGs.Enrichment <- rbind(DEGs.Enrichment,
                          data.frame(day = rep(names(DEGs)[i], times = length(T.MCLNum)),
                                     MCLNum = formatC(T.MCLNum, width = 3, flag = "0"),
                                     pvalue = T.pvalue,
                                     N = rep(N, times = length(T.MCLNum)),
                                     M = M,
                                     n = n,
                                     x = x,
                                     Intersection_AGI = paste(T.AGI, collapse = "|"), 
                                     stringsAsFactors = F)
  )
  print(length(DEGs)-i)
  i <- i+1
}
#save data----
saveRDS(object = DEGs.Enrichment, file = "~/nishioka/RDS/DEGs_Enrichment.rds")
write.table(x = DEGs.Enrichment, file = "~/nishioka/Table/DEGs_Enrichment.txt", sep = "\t", quote = F, row.names = F)
after <- proc.time()
print(after - before)#533.8 sec
#Input data----
DEGs.Enrichment <- readRDS(file = "~/nishioka/RDS/DEGs_Enrichment.rds")
MasterTable <- readRDS(file = "~/nishioka/RDS/MasterTable.rds")
#processing data----
df <- DEGs.Enrichment %>% select(-Intersection_AGI)
df <- df[as.numeric(df$MCLNum) < 70, ]
df$day <- factor(x = df$day,
                 levels = unique(df$day)[4:1],
                 ordered = TRUE)
#Divide p-value into five----
df$pvalue[df$pvalue > 5e-2] <- 1
df$pvalue[df$pvalue <= 5e-2 & df$pvalue > 5e-5] <- 5/6
df$pvalue[df$pvalue <= 5e-5 & df$pvalue > 5e-10] <- 4/6
df$pvalue[df$pvalue <= 5e-10 & df$pvalue > 5e-30] <- 3/6
df$pvalue[df$pvalue <= 5e-30 & df$pvalue > 5e-50] <- 2/6
df$pvalue[df$pvalue <= 5e-50 & df$pvalue > 5e-100] <- 1/6
df$pvalue[df$pvalue <= 5e-100] <- 0
#GSEV:plotting heatmap----
g <- ggplot(data = df, aes(x = MCLNum, y = day, fill = pvalue))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient(high="white",low="red", na.value = "white")
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g <- g + theme(axis.ticks.y = element_blank())
g <- g + theme(legend.position = "top")
g <- g + theme(legend.position = 'none')
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + theme(legend.title = element_blank())
g <- g + theme(legend.text = element_text(size=24))
ggsave(filename = "~/nishioka/Image/GSEV/DEGs_GSEV.png", plot = g, width = 8, height = 1)