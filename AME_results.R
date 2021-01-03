#package----
library(tidyverse)
#input data-----
#FDR001----
NodeTable <- read.table("~/IMO/base/DEG_FDR0.01/200916DEG_FDR0.01_PCC_FDR0.001_inflation2.5_iteration9_formastertable.csv", header=T, sep=",", stringsAsFactors = F)
T.filename <- list.files("~/IMO/Table/AME_upstream500/DEG_FDR0.01/Motif/Results/", full.names = T)
filename <- c()
Motif_ID <- c()
obnames <- c()
consensus <- c()
motif_alt_ID <- c()
T.MCLNum <- NodeTable$X__mclCluster %>% unique()
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
i <- 1
for(i in i:length(T.MCLNum)){
  filename <- T.filename[grep(paste0("MCLNum", T.MCLNum[i], "_"), T.filename)]
  title <- paste0(filename, "/ame.tsv")
  motif_ame_results <- try(read.table(title, sep = "\t", header = T, stringsAsFactors = F), silent = TRUE)
  if(class(motif_ame_results) != "try-error"){
    Motif_ID <- c(Motif_ID, motif_ame_results$motif_ID)
    consensus <- c(consensus, motif_ame_results$consensus)
    motif_alt_ID <- c(motif_alt_ID, motif_ame_results$motif_alt_ID)
    obnames <- c(obnames, rep(paste0("MCLNum", T.MCLNum[i]), times = length(motif_ame_results$motif_ID)))
  }
  i <- i+1
}

AME_results <- data.frame(MCLNum = obnames,
                          MotifID = Motif_ID,
                          consensu = consensus,
                          motif_alt_ID = motif_alt_ID,
                          stringsAsFactors = F
)
#AGRIS(https://agris-knowledgebase.org/AtcisDB/bindingsites.html)
#WRKY:TTGAC
#bHLH(G-box promoter motif):CACGTG
#ERF(GCC-box promoter motif):GCCGCC
#MYB binding site promoter:(A/C)ACC(A/T)A(A/C)C
#reference:Regulating the Regulators: The Control of Transcription Factors in Plant Defense Signaling
#bZIP:ACGT
#NAC:CGT[G/T], ACG
#RAV:CAACA
#Dof:AAAG
#SBP:GTAC
MotifList <- c("TTGAC", "GTCAA", "CACGTG", "GCCGCC", "GGCGGC", "ACGT", "CGT", "ACG", "CAACA", "TGTTG", "AAAG", "CTTT", "GTAC")
names(MotifList) <- c("WRKY", "WRKY", "bHLH", "ERF", "ERF", "bZIP", "NAC", "NAC", "RAV", "RAV", "Dof", "Dof", "SBP")
CoreMotif <- rep(0, times = nrow(AME_results))
TF <- rep(0, times = nrow(AME_results))
i <- 1
for(i in i:length(MotifList)){
  if(length(CoreMotif[grep(MotifList[i], AME_results$consensu)]) != 0){
    n <- 1
    for(n in n:length(CoreMotif[grep(MotifList[i], AME_results$consensu)])){
      if(!is.na(CoreMotif[grep(MotifList[i], AME_results$consensu)][n])){
        if(CoreMotif[grep(MotifList[i], AME_results$consensu)][n] ==  0){
          CoreMotif[grep(MotifList[i], AME_results$consensu)][n] <- MotifList[i]
          TF[grep(MotifList[i], AME_results$consensu)][n] <- names(MotifList)[i]
        }else{
          CoreMotif[grep(MotifList[i], AME_results$consensu)][n] <- paste0(CoreMotif[grep(MotifList[i], AME_results$consensu)][n], "|", MotifList[i])
          TF[grep(MotifList[i], AME_results$consensu)][n] <- paste0(TF[grep(MotifList[i], AME_results$consensu)][n], "|", names(MotifList)[i])
        }
      }
      n <- n+1
    }
  }
  i <- i+1
}
AME_results <-AME_results %>% mutate(CoreMotif = CoreMotif,
                                     TF = TF)
#save----
saveRDS(AME_results, file = "~/IMO/RDS/AME_results_FDR001_inflations2.5.rds")
write.table(AME_results, file = "~/IMO/Table/AME_results_FDR001_inflations2.5.txt", append=F, quote = F, sep = "\t", row.names = F)
#FDR005----
#input data-----
T.inflations <- c("inflations2.5", "inflations3.5", "inflations4.5")
FDR005.path <- list.files(path = "~/IMO/base/DEG_FDR0.05/")
h <- 1
for (h in h:length(T.inflations)) {
  T.filename <- list.files(paste0("~/IMO/Table/AME_upstream500/DEG_FDR0.05/", T.inflations[h], "/Motif/Results/"), full.names = T)
  T.path <- list.files(path = paste0("~/IMO/base/DEG_FDR0.05/", FDR005.path[1]))
  T.path <- T.path[grep("csv", T.path)]
  NodeTable <- read.table(paste0("~/IMO/base/DEG_FDR0.05/", FDR005.path[1], "/", T.path), header=T, sep=",", stringsAsFactors = F)
  filename <- c()
  Motif_ID <- c()
  obnames <- c()
  consensus <- c()
  motif_alt_ID <- c()
  T.MCLNum <- NodeTable$X__mclCluster %>% unique()
  T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
  i <- 1
  for(i in i:length(T.MCLNum)){
    filename <- T.filename[grep(paste0("MCLNum", T.MCLNum[i], "_"), T.filename)]
    title <- paste0(filename, "/ame.tsv")
    motif_ame_results <- try(read.table(title, sep = "\t", header = T, stringsAsFactors = F), silent = TRUE)
    if(class(motif_ame_results) != "try-error"){
      Motif_ID <- c(Motif_ID, motif_ame_results$motif_ID)
      consensus <- c(consensus, motif_ame_results$consensus)
      motif_alt_ID <- c(motif_alt_ID, motif_ame_results$motif_alt_ID)
      obnames <- c(obnames, rep(paste0("MCLNum", T.MCLNum[i]), times = length(motif_ame_results$motif_ID)))
    }
    i <- i+1
  }
  AME_results <- data.frame(MCLNum = obnames,
                            MotifID = Motif_ID,
                            consensu = consensus,
                            motif_alt_ID = motif_alt_ID,
                            stringsAsFactors = F
  )
  #AGRIS(https://agris-knowledgebase.org/AtcisDB/bindingsites.html)
  #WRKY:TTGAC
  #bHLH(G-box promoter motif):CACGTG
  #ERF(GCC-box promoter motif):GCCGCC
  #MYB binding site promoter:(A/C)ACC(A/T)A(A/C)C
  #reference:Regulating the Regulators: The Control of Transcription Factors in Plant Defense Signaling
  #bZIP:ACGT
  #NAC:CGT[G/T], ACG
  #RAV:CAACA
  #Dof:AAAG
  #SBP:GTAC
  MotifList <- c("TTGAC", "GTCAA", "CACGTG", "GCCGCC", "GGCGGC", "ACGT", "CGT", "ACG", "CAACA", "TGTTG", "AAAG", "CTTT", "GTAC")
  names(MotifList) <- c("WRKY", "WRKY", "bHLH", "ERF", "ERF", "bZIP", "NAC", "NAC", "RAV", "RAV", "Dof", "Dof", "SBP")
  CoreMotif <- rep(0, times = nrow(AME_results))
  TF <- rep(0, times = nrow(AME_results))
  i <- 1
  for(i in i:length(MotifList)){
    if(length(CoreMotif[grep(MotifList[i], AME_results$consensu)]) != 0){
      n <- 1
      for(n in n:length(CoreMotif[grep(MotifList[i], AME_results$consensu)])){
        if(!is.na(CoreMotif[grep(MotifList[i], AME_results$consensu)][n])){
          if(CoreMotif[grep(MotifList[i], AME_results$consensu)][n] ==  0){
            CoreMotif[grep(MotifList[i], AME_results$consensu)][n] <- MotifList[i]
            TF[grep(MotifList[i], AME_results$consensu)][n] <- names(MotifList)[i]
          }else{
            CoreMotif[grep(MotifList[i], AME_results$consensu)][n] <- paste0(CoreMotif[grep(MotifList[i], AME_results$consensu)][n], "|", MotifList[i])
            TF[grep(MotifList[i], AME_results$consensu)][n] <- paste0(TF[grep(MotifList[i], AME_results$consensu)][n], "|", names(MotifList)[i])
          }
        }
        n <- n+1
      }
    }
    i <- i+1
  }
  AME_results <-AME_results %>% mutate(CoreMotif = CoreMotif,
                                       TF = TF)
  #save----
  saveRDS(AME_results, file = paste0("~/IMO/RDS/AME_results_FDR005_", T.inflations[h], ".rds"))
  write.table(AME_results, file = paste0("~/IMO/Table/AME_results_FDR005_", T.inflations[h], ".txt"), append=F, quote = F, sep = "\t", row.names = F)
  h <- h+1
}
