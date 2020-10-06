#package---------------------------------------------------------
library(tidyverse)
#input data------------------------------------------------------
upstream_500 <- read.table(file = "~/IMO/PublishedData/TAIR10_upstream_500_20101028.txt", fill = T, sep = ",", stringsAsFactors = F)
upstream_500 <- as.character(unlist(upstream_500))
#processing data-------------------------------------------------
temp <- grep("chr", upstream_500)
T.AGI <- c()
i <- 1
total <- length(temp)
for(i in i:total){
  T.AGI <- c(T.AGI, substr(upstream_500[temp[i]], 2, 10))
  print(i)
  i <- i+1
}
#ÅŒã‚Ìˆêü‚¾‚¯©“®‰»‚Å‚«‚È‚©‚Á‚½B
allsequence <- list()
presequence <- c()
sequence <- c()
i <- 1
total <- length(temp)
for(i in i:c(total-1)){
  presequence <- upstream_500[c(temp[i]+1):c(temp[i+1]-1)]
  n <- 1
  for(n in n:length(presequence)){
    sequence <- paste0(sequence, presequence[n])
    n <- n+1
  }
  allsequence <- c(allsequence, list(sequence))
  sequence <- c()
  print(i)
  i <- i+1
}
#c‚è‚ÌÅŒã‚Ìˆêü‚ğ’Ç‰Á
presequence <- upstream_500[c(temp[total]+1):length(upstream_500)]
sequence <- paste0(sequence, presequence[n])
allsequence <- c(allsequence, list(sequence))
#data.frame
up_500bp <- data.frame(AGI = T.AGI,
                       sequence = unlist(allsequence),
                       stringsAsFactors = F
)
saveRDS(object = up_500bp, file = "~/IMO/RDS/TAIR10_upstream_500_20101028.rds")