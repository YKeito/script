#package----
library(tidyverse)
#input data-----
T.filename <- list.files("~/IMO/RDS/", full.names = T)
T.filename <- T.filename[grep("AME_results_FDR00", T.filename)]
T.FDR <- T.filename %>% str_split(pattern = "_", simplify = T)
T.inflations <- T.FDR[, 4] %>% str_sub(end = nchar(T.FDR[, 4])-4)
T.FDR <- paste0(T.FDR[, 3], "_", T.inflations)

T.inflations <- c("inflations2.5", "inflations3.5", "inflations4.5")
MotifList <- c("TTGAC", "GTCAA", "CACGTG", "GCCGCC", "GGCGGC", "ACGT", "CGT", "ACG", "CAACA", "TGTTG", "AAAG", "CTTT", "GTAC")
names(MotifList) <- c("WRKY", "WRKY", "bHLH", "ERF", "ERF", "bZIP", "NAC", "NAC", "RAV", "RAV", "Dof", "Dof", "SBP")

i <- 1
for (i in i:length(T.filename)) {
  AME_results <- readRDS(file = paste0(T.filename[i]))
  AME_results$TF[AME_results$TF == 0] <- "others"
  TFname <- c("others", unique(names(MotifList)))
  allMCLNum <- list()
  j <- 1
  for(j in j:length(TFname)){
    T_MotifID <- AME_results$MotifID[grep(TFname[j], AME_results$TF)]
    T.MCLNum <- c()
    k <- 1
    for(k in k:length(T_MotifID)){
      T.MCLNum <- c(T.MCLNum, AME_results$MCLNum[T_MotifID[k] == AME_results$MotifID])
      k <- k+1
    }
    allMCLNum <- c(allMCLNum, list(T.MCLNum))
    j <- j+1
  }
  names(allMCLNum) <- TFname
  UniMCLNum <- unique(AME_results$MCLNum)
  count <- list()
  j <- 1
  for(j in j:length(UniMCLNum)){
    total <- c()
    k <- 1
    for(k in k:length(allMCLNum)){
      total <- c(total, sum(allMCLNum[[k]] == UniMCLNum[j], na.rm = T)/sum(unlist(allMCLNum) == UniMCLNum[j], na.rm = T))
      k <- k+1
    }
    names(total) <- TFname
    count <- c(count, list(total))
    j <- j+1
  }
  names(count) <- UniMCLNum
  temp <- str_sub(UniMCLNum, start = 7, end = 9)
  T_data <- data.frame(MCLNum = rep(paste0(str_sub("0000", start = 1, end  = nchar("000")-nchar(temp)), temp), each = length(TFname)),
                       value = as.vector(unlist(count)),
                       TF = rep(paste0("0", length(TFname):1, TFname), times = length(UniMCLNum)),
                       stringsAsFactors = F
  )
  g <- ggplot(T_data, aes(x = TF, y = MCLNum, fill = value))
  g <- g + geom_tile(color = "black", size = 0.1)
  g <- g + scale_fill_gradient2(high = "red")
  g <- g+theme_dark()
  g <- g + theme_linedraw()
  g <- g + coord_flip()
  g <- g +@theme(axis.text=element_text(size=15))
  g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 0))
  title <- paste0("~/IMO/Image/AME/", T.FDR[i], "_enrichment_heatmap.png")
  
  ggsave(title, g, width = 8, height = 5)  
  i <- i+1
}