####"C:/Users/Matsnaga-21/Desktop/R script (nishioka)/"####



data<- read.table("200827coexptable_possitive_0.05(FDR0.01).txt-inflation2.5_iteration100-clustered default node.csv", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
data[is.na(data)]<-0

n <- 1
total <- 50 #max(data$X__mclCluster)
count <- c()      
for(n in n:total){
  x <- length(which(data$X__mclCluster == n))
  count[n]<-x
  n <- n+1
}

tiff("200830_1-50_q0.0001_inflation5.5.tiff")
barplot(count, xlab = "MCLNum", ylab = "count", main = "1-50_q0.0001_inflation5.5")
dev.off()