#package----
library(tidyverse)
#input data----
#各MCLNumのAGIの配列をクラスター単位で引っこ抜く----
up_500bp <- readRDS(file = "~/IMO/RDS/TAIR10_upstream_500_20101028.rds")
NodeTable <- read.table("~/IMO/base/DEG_FDR0.05/inflation2.5(defalt)/200916DEG_FDR0.05_PCC_FDR0.0001_inflation2.5_iteration8_formastertable.csv", header=T, sep=",", stringsAsFactors = F)
T.MCLNum <- NodeTable$X__mclCluster %>% unique()
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
i <- 1
for(i in i:length(T.MCLNum)){
  T.AGI <- NodeTable %>% filter(X__mclCluster == T.MCLNum[i]) %>% select(name) %>% unlist(use.names = F)
  T.data <- up_500bp[match(T.AGI, up_500bp$AGI), ]
  T.Control <- sample(up_500bp$AGI, length(T.AGI))
  T.ControlData <- up_500bp[match(T.Control, up_500bp$AGI), ]
  fasta <- c()
  allfasta <- c()
  cont_fasta <- c()
  cont_allfasta <- c()
  total_o <- nrow(T.data)
  o <- 1
  for(o in o:total_o){
    fasta <- rbind(str_sub(T.data[, "sequence"][o], start=1, end=80),
                   str_sub(T.data[, "sequence"][o], start=81, end=160),
                   str_sub(T.data[, "sequence"][o], start=161, end=240),
                   str_sub(T.data[, "sequence"][o], start=241, end=320),
                   str_sub(T.data[, "sequence"][o], start=321, end=400),
                   str_sub(T.data[, "sequence"][o], start=401, end=480),
                   str_sub(T.data[, "sequence"][o], start=481, end=500)
    )
    cont_fasta <- rbind(str_sub(T.ControlData[, "sequence"][o], start=1, end=80),
                        str_sub(T.ControlData[, "sequence"][o], start=81, end=160),
                        str_sub(T.ControlData[, "sequence"][o], start=161, end=240),
                        str_sub(T.ControlData[, "sequence"][o], start=241, end=320),
                        str_sub(T.ControlData[, "sequence"][o], start=321, end=400),
                        str_sub(T.ControlData[, "sequence"][o], start=401, end=480),
                        str_sub(T.ControlData[, "sequence"][o], start=481, end=500)
    )
    data_fastaAGI <- paste0(">", T.data[, "AGI"][o])
    control_fastaAGI <- paste0(">", T.ControlData[, "AGI"][o])
    allfasta <- c(allfasta, rbind(data_fastaAGI, fasta))
    cont_allfasta <- c(cont_allfasta, rbind(control_fastaAGI, cont_fasta))
    o <- o+1
  }
  target <- paste0("~/IMO/Table/DEG_FDR0.05/inflations2.5/Motif/MultiFasta/Target/", "MCLNum", T.MCLNum[i], "_upstream500.fasta")
  control <- paste0("~/IMO/Table/DEG_FDR0.05/inflations2.5/Motif/MultiFasta/Control/", "control_MCLNum", T.MCLNum[i], "_upstream500.fasta")
  write.table(allfasta, file = target, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(cont_allfasta, file = control, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  print(i)
  i <- i+1
}
#remove object----
rm(list = ls())