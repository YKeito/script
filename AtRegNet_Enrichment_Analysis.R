#test
#git
"~/bigdata/yasue/nishioka/script/AtRegNet_Enrichment_Analysis.R"
before <- proc.time()
#load package----
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
#Input data----
MasterTable <- readRDS(file = "~/bigdata/yasue/nishioka/RDS/MasterTable.rds")
TAIR10 <- readRDS("~/bigdata/yasue/TAIR10_ShortName.rds")
AtRegNet <- read.table("~/bigdata/yasue/AGRIS/AtRegNet/20190729AtRegNet_modified.txt", sep = "\t", quote = "", fill = T, stringsAsFactors = F, header = T)
AtRegNet$TFLocus <- toupper(AtRegNet$TFLocus)
AtRegNet$TargetLocus <- toupper(AtRegNet$TargetLocus)
AtRegNet <- AtRegNet[nchar(AtRegNet$TFLocus) == 9, ]
TF.AGI <- AtRegNet %>% distinct(TFLocus) %>% unlist(use.names = F)
AGRIS <- read.table(file = "~/bigdata/yasue/AGRIS/AGRIS_TFLIST/families_data.tbl", sep = "\t", header = F, quote = "", stringsAsFactors = F)
AGRIS <- data.frame(AGI = toupper(AGRIS$V2),
                    AGRIS = AGRIS$V1,
                    stringsAsFactors = F)
AGRIS <- AGRIS %>% distinct(AGI, .keep_all = T)
PLAZA <- read.table(file = "~/bigdata/yasue/TFList/Ath_TF_list.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
PLAZA <- data.frame(AGI = toupper(PLAZA$Gene_ID),
                    PLAZA = PLAZA$Family,
                    stringsAsFactors = F)
PLAZA <- PLAZA %>% distinct(AGI, .keep_all = T)
PlantTF <- read.table(file = "~/bigdata/yasue/TFList/PlnTFDB_tf_1429098561.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
PlantTF <- data.frame(AGI = toupper(str_sub(PlantTF$Protein.ID, end = 9)),
                      PlantTF = PlantTF$Family,
                      PlantTF_Category = PlantTF$Category,
                      stringsAsFactors = F)
PlantTF <- PlantTF %>% distinct(AGI, .keep_all = T)
#Enrichment analysis for TF Target genes----
N <- TAIR10 %>% nrow()
T.MCLNum <- MasterTable$MCLNum[!is.na(MasterTable$MCLNum)] %>% unique()
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
    T.data <- MasterTable %>% filter(MCLNum == T.MCLNum[j])
    n <- c(n, nrow(T.data))
    Cluster.AGI <- T.data$AGI
    T.AGI <- c(T.AGI, intersect(Target.AGI, Cluster.AGI))
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
                                     Intersection_AGI = paste(T.AGI, collapse = "|"), 
                                     stringsAsFactors = F)
  )
  print(length(TF.AGI)-i)
  i <- i+1
}
AtRegNet.Table <- AtRegNet.Table %>%
  left_join(y = TAIR10, by = "AGI") %>%
  left_join(y = AGRIS, by = "AGI") %>%
  left_join(y = PLAZA, by = "AGI") %>%
  left_join(y = PlantTF, by = "AGI")
#save data----
saveRDS(object = AtRegNet.Table, file = "~/bigdata/yasue/nishioka/RDS/AtRegNetTable.rds")
write.table(x = AtRegNet.Table, file = "~/bigdata/yasue/nishioka/Table/AtRegNetTable.txt", sep = "\t", quote = F, row.names = F)
write.table(x = AtRegNet.Table %>% filter(grepl("Dof", AGRIS), pvalue < 0.05, M >= 10), 
            file = "~/bigdata/yasue/nishioka/Table/DofTF_AtRegNetTable.txt", sep = "\t", quote = F, row.names = F)
after <- proc.time()
print(after - before)#533.8 sec