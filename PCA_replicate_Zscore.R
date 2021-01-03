#local R-studio version
#"~/nishioka/script/PCA_replicate_Zscore.R"
#install.packages("ggrepel")
#install.packages("tidyverse")
#load library----
library(ggplot2)
library(ggrepel)
library(stringr)
library(tidyverse)
library(factoextra)
#Input published data & formatting data----
Nishioka.data <- read.table(file = "~/nishioka/base/20200317_log2FC&FDR.txt", sep = "\t", header = T, stringsAsFactors = F)
Nishioka.ppm <- read.table(file = "~/nishioka/base/20200317_ExpTable.txt", sep = "\t", header = T, stringsAsFactors = F)
callus.data <- read.table(file = "~/nishioka/PublishedData/RNA-Seq_callus.txt", sep = "\t", header = T, stringsAsFactors = F)
#detecting DEGs----
DEGs <- list()
T.names <- c("2d", "4d", "7d", "10d")
df <- c()
i <- 1
for(i in i:length(T.names)){
  Nishioka.DEGs <- Nishioka.data %>% select(contains(T.names[i]))
  FDR.data <- Nishioka.DEGs %>% select(ends_with("qvalue"))
  FDR.logic <- FDR.data < 0.05
  FDR.logic <- FDR.logic %>% apply(MARGIN = 1, FUN = sum) == 2
  DEGs <- c(DEGs, list(Nishioka.data$reference_id[FDR.logic]))
  log2FC.data <- Nishioka.DEGs %>% select(ends_with("log2FC"))
  Up.logic <- log2FC.data > 0 & FDR.logic
  Down.logic <- log2FC.data < 0 & FDR.logic
  df <- rbind(df,
              data.frame(Sample = T.names[i],
                         Group = c("02Down", "01Up"),
                         value = c(sum(Down.logic), sum(Up.logic)))
  )
  i <- i+1
}
names(DEGs) <- T.names
DEGs.AGI <- DEGs %>% unlist() %>% unique()
#VennDiagram----
data <- DEGs
Venn <- venn(data)
names(attr(Venn,"intersections")) <- paste0(names(attr(Venn,"intersections")), ",")
T.Venn <- data.frame(GeneList = str_split(unlist(attr(Venn,"intersections"), use.names = F), pattern = ",", simplify = T),
                     conditions = str_split(names(unlist(attr(Venn,"intersections"))), pattern = ",", simplify = T)[, 1],
                     stringsAsFactors = F)
saveRDS(object = T.Venn, file = "~/nishioka/RDS/Nishioka_VennDiagram.rds")
write.table(x = T.Venn, file = "~/nishioka/Table/Nishioka_VennDiagram.txt", sep = "\t", quote = F, row.names = F)
ggsave(filename = "~/nishioka/Image/DEGs_Graph/Nishioka_VennDiagram.png", plot = plot(Venn), width = 6, height = 4)
#DEGs;BarGraph----
df$Sample <- factor(x = df$Sample,
                    levels = df$Sample[1:8],
                    ordered = TRUE)
g <- ggplot(df, aes(x = Sample, y = value, fill = Group))
g <- g + geom_bar(stat = "identity")
g <- g + theme_bw()
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) #?????Á‚?
g <- g + theme(legend.position = 'none') #?}???Á‚?
g <- g + theme(axis.text = element_text(size = 18))
plot(g)
ggsave(filename = "~/nishioka/Image/DEGs_Graph/NumDEGs_vs0_BarGraph.png", plot = g, width = 8, height = 6)
#Convert ppm to Zscore----
#look up DEGs
DEGs.ppm <- Nishioka.ppm %>% filter(reference_id %in% DEGs.AGI)
#join published data----
T.data <- DEGs.ppm
total <- T.data %>% ncol()
Normalized.data <- data.frame(AGI = T.data$reference_id,
                              stringsAsFactors = F)
j <- 2
for (j in j:total) {
  Normalized.data <- cbind(Normalized.data, 
                           T.data %>% select(j) %>% scale())
  j <- j+1
}
#PCA Replicate ver----
PCA.data <- Normalized.data %>% select(-AGI)
PCA.data <- PCA.data %>% na.omit()
PCA.data <- PCA.data %>% lapply(FUN = as.numeric) %>% data.frame()
PCA.data <- PCA.data %>% t()
PCA <- prcomp(PCA.data, scale = TRUE)
T.component <- summary(PCA)
T.component <- T.component$importance[2, 1:2]
#PCA;point plot----
df <- data.frame(PC1 = PCA$x[, 1],
                 PC2 = PCA$x[, 2],
                 Group = str_split(rownames(PCA$x), pattern = "_log", simplify = T)[, 1],
                 row.names = NULL
)
df <- df %>% mutate(Sample = str_split(df$Group, pattern = "_", simplify = T)[, 1])
df <- df %>% mutate(Time = str_split(df$Group, pattern = "_", simplify = T)[, 2])
g <- ggplot(data = df, aes(x=PC1, y=PC2, color = Time, label = Group, shape = Sample))
g <- g + theme_bw()
g <- g + geom_point(size = 5)
g <- g + xlab(paste0("PC1 (", formatC(T.component[1]*100, width = 3), "%)"))
g <- g + ylab(paste0("PC2 (", formatC(T.component[2]*100, width = 3), "%)"))
g <- g + theme(axis.text = element_text(size=18))
g <- g + theme(axis.title = element_text(size=18))
ggsave(filename = "~/nishioka/Image/PCA/NishiokaData_replicate_PCA.png", plot = g, width = 9, height = 7)
#conver ppm to Average----
T.data <- DEGs.ppm %>% 
  left_join(y = callus.data, by = c("reference_id"="X"))
T.colnames <- T.data %>% colnames() %>% str_sub(end = -3) %>% unique()
T.colnames <- T.colnames[-1]
total <- T.colnames %>% length()
Ave.data <- data.frame(AGI = T.data$reference_id)
j <- 1
for (j in j:total) {
  Ave.data <- cbind(Ave.data, 
                    T.data %>% 
                      select(starts_with(T.colnames[j])) %>% 
                      apply(MARGIN = 1, FUN = mean)
  )
  j <- j+1
}
colnames(Ave.data)[2:ncol(Ave.data)] <- T.colnames
#conver Average to log2FC----
T.day <- c("day2", "day4", "day7", "day10")
T.variation <- c("Col_0", "rpt5a")
log2FC.data <- data.frame(AGI = Ave.data$AGI)
i <- 1
for (i in i:length(T.variation)) {
  j <- 1
  for (j in j:length(T.day)) {
    Numerator <- Ave.data %>% select(paste0(T.day[j], "_", T.variation[i]))+1
    denominator <- Ave.data %>% select(paste0("day0_", T.variation[i]))+1
    ratio <- Numerator/denominator
    log2FC.data <- cbind(log2FC.data,
                         ratio %>% log2()
    )
    j <- j+1
  }
  i <- i+1
}
Numerator <- Ave.data$Col_C14+1
denominator <- Ave.data$Col_C0+1
ratio <- Numerator/denominator
log2FC.data <- cbind(log2FC.data,
                     data.frame(callus = ratio %>% log2())
)
#conver log2FC to Zscore----
T.data <- log2FC.data
total <- T.data %>% ncol()
Normalized.log2FC <- data.frame(AGI = T.data$AGI,
                                stringsAsFactors = F)
j <- 2
for (j in j:total) {
  Normalized.log2FC <- cbind(Normalized.log2FC, 
                             T.data %>% select(j) %>% scale())
  j <- j+1
}
PCA.data <- Normalized.log2FC %>% select(-AGI)
PCA.data <- PCA.data %>% na.omit()
PCA.data <- PCA.data %>% lapply(FUN = as.numeric) %>% data.frame()
PCA.data <- PCA.data %>% t()
PCA <- prcomp(PCA.data, scale = TRUE)
T.component <- summary(PCA)
write.table(x = T.component$importance, file = "~/nishioka/Table/PCA_summary.txt", sep = "\t", quote = F, row.names = T)
write.table(x = PCA$x, file = "~/nishioka/Table/PCA_ScatterPlot_Coordinates.txt", sep = "\t", quote = F, row.names = T)
T.component <- T.component$importance[2, 1:2]
#PCA point plot;log2FC Zscore----
df <- data.frame(PC1 = PCA$x[, 1],
                 PC2 = PCA$x[, 2],
                 Group = rownames(PCA$x),
                 row.names = NULL
)
df <- df %>% mutate(Shape = c(rep("col_0", times = 4), rep("rpt5a", times = 4), "callus"))
df <- df %>% mutate(Color = c(rep(T.day, times = 2), "day14"))
g <- ggplot(data = df, aes(x=PC1, y=PC2, color = Color, label = Group, shape = Shape))
g <- g + theme_bw()
g <- g + geom_point(size = 10)
g <- g + xlab(paste0("PC1 (", formatC(T.component[1]*100, width = 3), "%)"))
g <- g + ylab(paste0("PC2 (", formatC(T.component[2]*100, width = 3), "%)"))
g <- g + theme(axis.text = element_text(size=20))
g <- g + theme(axis.title = element_text(size=20))
g <- g + theme(legend.text = element_text(size=12))
ggsave(filename = "~/nishioka/Image/PCA/Nishioka&Callus_PCA.png", plot = g, width = 9, height = 7)
#log2FC;heatmap----
T.data <- log2FC.data %>% select(-callus) %>% select(-contains("Col_0"))
df <- T.data %>% gather(key = "Sample", value = Exp, -AGI)
df$Exp <- as.numeric(df$Exp)
df$Exp[is.na(df$Exp)] <- 0
df$Exp[df$Exp > 2] <- 2
df$Exp[df$Exp < -2] <- -2
T.matrix <- T.data[, 2:ncol(T.data)]
rownames(T.matrix) <- T.data$AGI
res <- hclust(dist(T.matrix), method = "ward.D2")
df$AGI <- factor(x = df$AGI,
                 levels = res$labels[res$order],
                 ordered = TRUE)
df$Sample <- factor(x = df$Sample,
                    levels = unique(df$Sample),
                    ordered = TRUE)
g <- ggplot(data = df, aes(x = Sample, y = AGI, fill = Exp))
g <- g + theme_bw()
g <- g + geom_tile()
g <- g + scale_fill_gradient2(limits = c(-2, 2), low = "cyan", high = "yellow", mid = "black", na.value = "white")
g <- g + ylab("") + xlab("")
g <- g + theme(axis.title.y = element_text(size = 14))
g <- g + theme(strip.text.x = element_text(size = 12))
g <- g + theme(axis.title.x = element_blank())
g <- g + theme(axis.ticks.y = element_blank())
g <- g + theme(legend.title = element_blank())
g <- g + theme(legend.text = element_text(size=12))
g <- g + theme(axis.text.x = element_text(size=24))
g <- g + theme(axis.text.y = element_blank())
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "~/nishioka/Image/Heatmap/log2FC_Heatmap.png", plot = g, width = 6, height = 16)
Heatmap.matrix <- data.frame(AGI = res$labels[res$order],
                             stringsAsFactors = F)
Heatmap.matrix <- Heatmap.matrix %>% 
  left_join(y = log2FC.data, by = "AGI") %>% 
  left_join(y = Nishioka.data %>% select(-contains("logFC")), by = c("AGI" = "reference_id")) 
write.table(x = Heatmap.matrix,
            file = "~/nishioka/Table/Nishioka_DEGsExp_HCL.txt", 
            sep = "\t", quote = F, row.names = F)
#Normalized log2FC(Zscore);heatmap----
T.data <- Normalized.log2FC
df <- T.data %>% gather(key = "Sample", value = Exp, -AGI)
df$Exp <- as.numeric(df$Exp)
df$Exp[is.na(df$Exp)] <- 0
df$Exp[df$Exp > 2] <- 2
df$Exp[df$Exp < -2] <- -2
T.matrix <- T.data[, 2:ncol(T.data)]
rownames(T.matrix) <- T.data$AGI
res <- hclust(dist(T.matrix), method = "ward.D2")
df$AGI <- factor(x = df$AGI,
                 levels = res$labels[res$order],
                 ordered = TRUE)
res.MCL <- hclust(dist(T.matrix %>% t()), method = "ward.D2")
df$Sample <- factor(x = df$Sample,
                    levels = res.MCL$labels[res.MCL$order],
                    ordered = TRUE)
g <- ggplot(data = df, aes(x = Sample, y = AGI, fill = Exp))
g <- g + theme_bw()
g <- g + geom_tile()
g <- g + scale_fill_gradient2(limits = c(-2, 2), low = "cyan", high = "yellow", mid = "black", na.value = "white")
g <- g + ylab("") + xlab("")
g <- g + theme(axis.title.y = element_text(size = 14))
g <- g + theme(strip.text.x = element_text(size = 12))
g <- g + theme(axis.title.x = element_blank())
g <- g + theme(axis.ticks.y = element_blank())
g <- g + theme(legend.title = element_blank())
g <- g + theme(legend.text = element_text(size=12))
g <- g + theme(axis.text.x = element_text(size=24))
g <- g + theme(axis.text.y = element_blank())
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "~/nishioka/Image/Heatmap/Normalized_log2FC_Zscore_Heatmap.png", plot = g, width = 12, height = 18)
#plotti