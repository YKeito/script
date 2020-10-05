#"C:/Users/Matsnaga-21/Desktop/R script (nishioka)/""
#input data-----------------------------------------------------------
allgenes <- read.table("190407 ver6/190407pvalueqvalue/190407Rfoldchange&p&q.txt", fill = T, quote = "", sep = "\t", header = T, stringsAsFactors = F)
annotation <- allgenes$annotation
names(annotation) <- allgenes$reference_id
allgenes <- allgenes$reference_id
NodeTable <- read.table("190409q0.0001inflation4.5iteration8.txt", header=T, sep="\t", stringsAsFactors = F)
TFList <- read.table("190222TF_list.txt", fill = T, quote = "", sep="\t", stringsAsFactors = F)
allRNASeq <- read.table("190407FDR0.05.txt", fill = T, quote = "", header=T, sep="\t", stringsAsFactors = F)
genesymbol <- read.table("summary.csv", fill = T, quote = "", sep = ",", header = T, stringsAsFactors = F) #昇順にしておくこと
genesymbol$input <- substr(genesymbol$input, start = 1, stop = c(nchar(genesymbol$input)-1)) #半角スペース入っている時だけ必要

#AGI--------------------------------------------------------------
#遺伝子リストを用意
colnames(TFList) <- c("TF_family", "AGI", "TF_name")
TFList$AGI <- toupper(TFList$AGI)

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