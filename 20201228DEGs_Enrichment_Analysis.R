#load package----
library(tidyverse)
#Input data----
TAIR10 <- readRDS(file = "~/IMO/RDS/TAIR10_ShortName.rds")
MCL.data <- read.table(file = "~/IMO/base/200916DEG_FDR0.05_PCC_FDR0.0001_inflation4.5_iteration10_formastertable.csv", sep = ",", header = T, fill = T)
MCL.data <- MCL.data %>% select(X__mclCluster, name)
DEGs.Venn <- read.table(file = "~/IMO/Table/20201217Nishioka_VennDiagram.txt", sep = "\t", header = T, quote = "")
#detecting DEGs----
T.DEGs <- list()
T.names <- c("02d", "04d", "07d", "10d")
i <- 1
for(i in i:length(T.names)){
  T.DEGs <- c(T.DEGs, list(DEGs.Venn %>% 
                             filter(conditions %in% T.names[i]) %>% 
                             select(GeneList) %>% 
                             unlist(use.names = F)
                           )
              )
}
names(T.DEGs) <- T.names
#Enrichment analysis for TF Target genes----
N <- TAIR10 %>% nrow()
T.MCLNum <- MCL.data$X__mclCluster %>% unique()
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
DEGs.Enrichment <- c()
i <- 1
for(i in i:length(T.DEGs)){
  Target.AGI <- T.DEGs[[i]]
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
  DEGs.Enrichment <- rbind(DEGs.Enrichment,
                           data.frame(day = rep(names(T.DEGs)[i], times = length(T.MCLNum)),
                                      MCLNum = formatC(T.MCLNum, width = 3, flag = "0"),
                                      pvalue = T.pvalue,
                                      N = rep(N, times = length(T.MCLNum)),
                                      M = M,
                                      n = n,
                                      x = x,
                                      Intersection_AGI = T.AGI, 
                                      stringsAsFactors = F)
  )
  print(length(T.DEGs)-i)
  i <- i+1
}
#save data----
saveRDS(object = DEGs.Enrichment, file = "~/IMO/RDS/20201228DEGs_Enrichment.rds")
write.table(x = DEGs.Enrichment, file = "~/IMO/Table/20201228DEGs_Enrichment.txt", sep = "\t", quote = F, row.names = F)