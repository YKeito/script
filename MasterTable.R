#package----
library(tidyverse)
#input data-----------------------------------------------------------
Nishioka.data <- read.table(file = "~/IMO/base/20200317_log2FC&FDR.txt", sep = "\t", header = T, stringsAsFactors = F)
NodeTable <- read.table("~/IMO/base/DEG_FDR0.01/200916DEG_FDR0.01_PCC_FDR0.001_inflation2.5_iteration9_formastertable.csv", header=T, sep=",", stringsAsFactors = F)
AGRIS <- read.table("~/IMO/PublishedData/AtTFDB/families_data.tbl", fill = T, quote = "", sep="\t", stringsAsFactors = F) %>% select(1, 2)
colnames(AGRIS) <- c("TF_family", "AGI")
AGRIS <- AGRIS %>% distinct(AGI, .keep_all = T)
TAIR10 <- read.table(file = "~/IMO/PublishedData/Arabidopsis_thaliana.TAIR10.48.gtf", sep = "\t", fill = T)
TAIR10 <- TAIR10 %>% filter(V3 == "gene") %>% select(V9)
TAIR10 <- TAIR10$V9 %>% str_split(pattern = ";", simplify = T)
GeneID <- TAIR10[, 1] %>% str_sub(start = 9)
ShortName <- TAIR10[, 2] %>% str_sub(start = 12)
ShortName[ShortName == "e araport11"] <- GeneID[ShortName == "e araport11"]
TAIR10.ShortName <- data.frame(AGI = GeneID,
                               ShortName = ShortName,
                               stringsAsFactors = F)
saveRDS(TAIR10.ShortName, "~/IMO/RDS/TAIR10_ShortName.rds")
#MasterTable作成本番-------------------------------------------------------------------
time <- c("X0_2d", "X0_4d", "X0_7d", "X0_10d")
T.qvalue <- paste0(time, "_qvalue")
allsample <- c()
obnames <- c()
n <- 1
for(n in n:length(time)){
  temp <- data.frame(AGI = rownames(allRNASeq[allRNASeq[, T.qvalue[n]] < 0.05, ]), stringsAsFactors = F)
  T_data <- rep("No", times = length(allgenes))
  names(T_data) <- allgenes
  T_data[match(temp[, "AGI"], names(T_data))] <- "Yes"
  allsample <- cbind(allsample, T_data)
  obnames <- c(obnames, time[n])
  n <- n+1
}
colnames(allsample) <- paste0(obnames, "_FDR005")

#AGI, TF-----------------------------------------------------------------------------------------
T.TF <- rep("No", times = length(allgenes))
names(T.TF) <- allgenes
T.AGI <- intersect(TFList$AGI, names(T.TF))
T.TF[match(T.AGI, names(T.TF))] <- TFList$TF_family[match(T.AGI, TFList$AGI)]

#MCL--------------------------------------------------------------------------------------------
MCLNum <- rep(NA, times = length(allgenes))
names(MCLNum) <- allgenes
MCLNum[match(NodeTable$name, names(MCLNum))] <- NodeTable$X__mclCluster

#Degree-------------------------------------------------------------------
Degree <- rep(0, times = length(allgenes))
names(Degree) <- allgenes
Degree[match(NodeTable$name, names(Degree))] <- NodeTable$Degree


#genesymbol
G.genesymbol <- rep("No", times = length(allgenes))
names(G.genesymbol) <- allgenes
G.AGI <- intersect(genesymbol$input, names(G.genesymbol))
G.genesymbol[match(G.AGI, names(G.genesymbol))] <- genesymbol$symbol[match(G.AGI, genesymbol$input)]
G.genesymbol[match(genesymbol$input, names(G.genesymbol))] <- genesymbol$symbol


#BetweennessCentrality-------------------------------------------------------------------
BC <- rep(0, times = length(allgenes))
names(BC) <- allgenes
BC[match(NodeTable$name, names(BC))] <- NodeTable$BetweennessCentrality
#MasterTable-------------------------------------------------------------------
MasterTable <- data.frame(AGI = allgenes,
                          TF = T.TF,
                          MCLNum = MCLNum,
                          Degree = Degree,
                          BetweennessCentrality = BC,
                          allsample,
                          genesymbol = genesymbol$symbol,
                          annotation = annotation,
                          stringsAsFactors = F
)
rownames(MasterTable) <- c()
#check
head(na.omit(MasterTable[MasterTable$MCLNum == 44, ]))
#save-------------------------------------------------------------------
save(MasterTable, file = "190409fdr0.05MCLmastertable.RData")
write.table(MasterTable, file = "190409fdr0.05MasterTable.txt", append=F, quote = F, sep = "\t", row.names = F)