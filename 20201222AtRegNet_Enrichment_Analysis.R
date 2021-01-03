before <- proc.time()
#load library----
library(tidyverse)
#Input data----
TAIR10 <- readRDS(file = "~/IMO/RDS/TAIR10_ShortName.rds")
MCL.data <- read.table(file = "~/IMO/base/200916DEG_FDR0.05_PCC_FDR0.0001_inflation4.5_iteration10_formastertable.csv", sep = ",", header = T, fill = T)
MCL.data <- MCL.data %>% select(X__mclCluster, name)
AtRegNet <- read.table("~/IMO/PublishedData/20200418AtRegNet_modified.txt", sep = "\t", quote = "", fill = T, stringsAsFactors = F, header = T)
AtRegNet$TFLocus <- toupper(AtRegNet$TFLocus)
AtRegNet$TargetLocus <- toupper(AtRegNet$TargetLocus)
AtRegNet <- AtRegNet[nchar(AtRegNet$TFLocus) == 9, ]
TF.AGI <- AtRegNet %>% distinct(TFLocus) %>% unlist(use.names = F)
#Enrichment analysis for TF Target genes----
N <- TAIR10 %>% nrow()
T.MCLNum <- MCL.data$X__mclCluster %>% unique()
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
AtRegNet.Table <- c()
i <- 1
for(i in i:length(TF.AGI)){
  Target.AGI <- AtRegNet %>% filter(TFLocus == TF.AGI[i]) %>% select(TargetLocus) %>% unlist(use.names = F) %>% unique()
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
  AtRegNet.Table <- rbind(AtRegNet.Table,
                          data.frame(AGI = rep(TF.AGI[i], times = length(T.MCLNum)),
                                     AtRegNet_TFFamily = rep(AtRegNet$TFName[match(TF.AGI[i], AtRegNet$TFLocus)], times = length(T.MCLNum)),
                                     MCLNum = formatC(T.MCLNum, width = 3, flag = "0"),
                                     pvalue = T.pvalue,
                                     N = rep(N, times = length(T.MCLNum)),
                                     M = M,
                                     n = n,
                                     x = x,
                                     Intersection_AGI = T.AGI, 
                                     stringsAsFactors = F)
  )
  print(length(TF.AGI)-i)
  i <- i+1
}
AtRegNet.Table <- AtRegNet.Table %>%
  left_join(y = TAIR10, by = "AGI")
#save data----
saveRDS(object = AtRegNet.Table, file = "~/IMO/RDS/20201222AtRegNetTable.rds")
write.table(x = AtRegNet.Table, file = "~/IMO/Table/20201222AtRegNetTable.txt", sep = "\t", quote = F, row.names = F)
after <- proc.time()
print(after - before)#533.8 sec