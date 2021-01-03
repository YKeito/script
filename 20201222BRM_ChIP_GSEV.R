#load library----
library(tidyverse)
#Input published data & formatting data----
DEGs.Ave <- read.table("~/IMO/base/20201217_Aveppm_DEGs.txt", sep = "\t", header = T, quote = "")
Nishioka.data <- read.table("~/IMO/base/20201217_log2FC.txt", sep = "\t", header = T, quote = "")
Nishioka.data <- Nishioka.data %>% filter(reference_id %in% DEGs.Ave$AGI)
Nishioka.data <- Nishioka.data %>% select(-6:-18)
TAIR10 <- readRDS(file = "~/IMO/RDS/TAIR10_ShortName.rds")
MCL.data <- read.table(file = "~/IMO/base/200916DEG_FDR0.05_PCC_FDR0.0001_inflation4.5_iteration10_formastertable.csv", sep = ",", header = T, fill = T)
MCL.data <- MCL.data %>% select(X__mclCluster, name)
BRM.ChIP <- read.table(file = "~/IMO/PublishedData/BRM_ChIP.txt", sep = "\t", header = T, quote = "", fill = T)
colnames(BRM.ChIP) <- "AGI"
#Enrichment analysis for TF Target genes----
N <- TAIR10 %>% nrow()
T.MCLNum <- MCL.data$X__mclCluster %>% unique()
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
Target.AGI <- BRM.ChIP$AGI
M <- Target.AGI %>% length()
T.pvalue <- c()
T.AGI <- c()
n <- c()
x <- c()
j <- 1
for(j in j:length(T.MCLNum)){
  T.data <- MCL.data %>% filter(X__mclCluster == T.MCLNum[j])
  n <- c(n, nrow(T.data))
  Cluster.AGI <- T.data$name
  T.AGI <- c(T.AGI, paste(intersect(Target.AGI, Cluster.AGI), collapse = "|"))
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
                       Intersection_AGI = T.AGI, 
                       stringsAsFactors = F
)
#save data----
saveRDS(object = BRM.GSEV, file = "~/IMO/RDS/20201222BRM_GSEV.rds")
write.table(x = BRM.GSEV, file = "~/IMO/Table/20201222BRM_GSEV.txt", sep = "\t", quote = F, row.names = F)
write.table(x = BRM.GSEV %>% filter(pvalue < 0.05), 
            file = "~/IMO/Table/20201222BRM_GSEV_pvalue005.txt", sep = "\t", quote = F, row.names = F)