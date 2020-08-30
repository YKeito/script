"~/nishioka/script/BRM_ChIP_GSEV.R"
#load library----
library(dplyr)
library(stringr)
#Input data----
TAIR10 <- read.table(file = "~/nishioka/PublishedData/TAIR10.txt", sep = "\t", header = T)
TAIR10$AGI <- TAIR10$AGI %>% str_sub(end = -3)
TAIR10 <- TAIR10 %>% distinct(AGI)
MasterTable <- readRDS(file = "~/nishioka/RDS/MasterTable.rds")
BRM.ChIP <- read.table(file = "~/nishioka/PublishedData/BRM_ChIP.txt", sep = "\t", header = T)
colnames(BRM.ChIP) <- "AGI"
#Enrichment analysis for TF Target genes----
N <- TAIR10 %>% nrow()
T.MCLNum <- MasterTable$MCLNum[!is.na(MasterTable$MCLNum)] %>% unique()
i <- 1
Target.AGI <- BRM.ChIP$AGI
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
BRM.GSEV <- data.frame(MCLNum = formatC(T.MCLNum, width = 3, flag = "0"),
                       pvalue = T.pvalue,
                       N = rep(N, times = length(T.MCLNum)),
                       M = M,
                       n = n,
                       x = x,
                       Intersection_AGI = paste(T.AGI, collapse = "|"), 
                       stringsAsFactors = F
                       )
#save data----
saveRDS(object = BRM.GSEV, file = "~/nishioka/RDS/BRM_GSEV.rds")
write.table(x = BRM.GSEV, file = "~/nishioka/Table/BRM_GSEV.txt", sep = "\t", quote = F, row.names = F)
write.table(x = BRM.GSEV %>% filter(pvalue < 0.05), 
            file = "~/nishioka/Table/BRM_GSEV_pvalue005.txt", sep = "\t", quote = F, row.names = F)
