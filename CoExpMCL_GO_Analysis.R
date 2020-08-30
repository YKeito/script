"~/bigdata/yasue/nishioka/script/CoExpMCL_GO_Analysis.R"
before <- proc.time()
#load package----
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
#Input data----
MasterTable <- readRDS(file = "~/bigdata/yasue/nishioka/RDS/MasterTable.rds")
TAIR10 <- readRDS("~/bigdata/yasue/TAIR10_ShortName.rds")
#TAIR10.GO <- read.table(file = "~/bigdata/yasue/GO_analysis/tair.gaf", skip = 24, sep = "\t", quote = "", stringsAsFactors = F, fill = T)
#TAIR10.GO <- read.table(file = "~/bigdata/yasue/GO_analysis/ATH_GO_GOSLIM.txt", skip = 4, sep = "\t", quote = "", stringsAsFactors = F, fill = T)
TAIR10.GO <- read.table(file = "/root/db/ATH_GO/ATH_GO_GOSLIM.txt", skip = 4, sep = "\t", quote = "", stringsAsFactors = F, fill = T)
colnames(TAIR10.GO) <- c("locus name",
                         "TAIR accession",
                         "object name",
                         "relationship",
                         "GO term",
                         "GO ID",
                         "TAIR Keyword ID",
                         "Aspect",
                         "GOslim term",
                         "Evidence",
                         "Evidence description",
                         "Evidence with",
                         "Reference",
                         "Annotator",
                         "Data annotated")
TAIR10.GO <- TAIR10.GO %>% select(`locus name`, `GO term`, `GO ID`, Aspect, `GOslim term`)
#GO analysis----
N <- TAIR10 %>% nrow()
T.MCLNum <- MasterTable$MCLNum[!is.na(MasterTable$MCLNum)] %>% unique()
T.MCLNum <- T.MCLNum[T.MCLNum <= 69]
T.Aspect <- c("P", "C", "F")
allGO.Table <- c()
i <- 1
for(i in i:length(T.Aspect)){
  Aspect.data <- TAIR10.GO %>% filter(Aspect == T.Aspect[i])
  UniGO.ID <- Aspect.data$`GO ID` %>% unique()
  GO.Table <- c()
  j <- 1
  for(j in j:length(T.MCLNum)){
    T.data <- MasterTable %>% filter(MCLNum == T.MCLNum[j])
    M <- nrow(T.data)
    Cluster.AGI <- T.data$AGI
    T.pvalue <- c()
    T.AGI <- c()
    x <- c()
    n <- c()
    GO.term <- c()
    GOslim.term <- c()
    k <- 1
    for(k in k:length(UniGO.ID)){
      UniGO.data <- Aspect.data %>% filter(`GO ID` == UniGO.ID[k]) %>% distinct(`locus name`, .keep_all = T)
      n <- c(n, UniGO.data %>% nrow())
      GO.term <- c(GO.term, paste(UniGO.data$`GO term` %>% unique(), collapse = "|"))
      GOslim.term <- c(GOslim.term, paste(UniGO.data$`GOslim term` %>% unique(), collapse = " | "))
      UniGO.AGI <- UniGO.data$`locus name`
      T.AGI <- c(T.AGI, paste(intersect(UniGO.AGI, Cluster.AGI), collapse = " | "))
      x <- c(x, length(intersect(UniGO.AGI, Cluster.AGI)))
      T.pvalue <- c(T.pvalue, phyper(x[k]-1, M, N-M, n[k], lower.tail = F))
      k <- k+1
    }
    GO.Table <- rbind(GO.Table,
                      data.frame(Aspect = rep(T.Aspect[i], times = length(UniGO.ID)),
                                 GOID = UniGO.ID,
                                 GOTerm = GO.term,
                                 GOslimTerm = GOslim.term,
                                 MCLNum = rep(formatC(T.MCLNum[j], width = 3, flag = "0"), times = length(UniGO.ID)),
                                 pvalue = T.pvalue,
                                 FDR = p.adjust(p = T.pvalue, method = "BH"),
                                 N = rep(N, times = length(UniGO.ID)),
                                 M = rep(M, times = length(UniGO.ID)),
                                 n = n,
                                 x = x,
                                 AGI = T.AGI,
                                 stringsAsFactors = F)
    )
    print(length(T.MCLNum)-j)
    j <- j+1
  }
  allGO.Table <- rbind(allGO.Table,
                       GO.Table
  )
  i <- i+1
}
saveRDS(object = allGO.Table, file = "~/bigdata/yasue/nishioka/RDS/GOTable.rds")
write.table(x = allGO.Table, file = "~/bigdata/yasue/nishioka/Table/GOTable.txt", sep = "\t", quote = F, row.names = F)
after <- proc.time()
print(after - before)#1588.94 sec, 26.47 min