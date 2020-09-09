####DEGs (FDRcut)####

####read.table####
data <- read.table("GRN_input/02_DEGs (FDRcut)/200817foldchange&pvalue&qvalue.txt", header=TRUE, row.names = 1, sep = "\t", quote = "")

####setting param####
param <- 0.05

####FDRcut####
subdata <- c()
subdata <- data[data$X0_2d_qvalue < param & data$X2d_qvalue < param | data$X0_4d_qvalue < param & data$X4d_qvalue < param | data$X0_7d_qvalue < param & data$X7d_qvalue < param | data$X0_10d_qvalue < param & data$X10d_qvalue < param, ]

####output####
write.table(subdata,"GRN_output/02_DEGs (FDRcut)/200817DEGs(FDR0.05).txt", append=F, quote = F, sep = "\t", row.names = T)



