#local R-studio version
#"~/nishioka/script/PCA&heatmap_Zscore.R"
#install.packages("ggrepel")
#install.packages("tidyverse")
#load library----
library(ggplot2)
library(ggrepel)
library(stringr)
library(tidyverse)
#Input published data & formatting data----
Nishioka.data <- read.table(file = "~/nishioka/base/20200317_ExpTable.txt", sep = "\t", header = T, stringsAsFactors = F)
callus.data <- read.table(file = "~/nishioka/PublishedData/RNA-Seq_callus.txt", sep = "\t", header = T, stringsAsFactors = F)
#Convert expression data to Zscore
T.colnames <- Nishioka.data %>% colnames() %>% str_sub(end = -3)
T.colnames <- T.colnames %>% unique()
T.colnames <- T.colnames[-1]
total <- T.colnames %>% length()
Normalized.Noshioka.data <- data.frame(AGI = Nishioka.data$reference_id,
                                       stringsAsFactors = F)
j <- 1
for (j in j:total) {
  Normalized.Noshioka.data <- cbind(Normalized.Noshioka.data, 
                                    Nishioka.data %>% 
                                      select(starts_with(T.colnames[j])) %>% 
                                      apply(MARGIN = 1, FUN = mean) %>% 
                                      scale()
                                    )
  j <- j+1
}
colnames(Normalized.Noshioka.data)[2:ncol(Normalized.Noshioka.data)] <- T.colnames
#callus data
T.colnames <- callus.data %>% colnames() %>% str_sub(end = -3)
T.colnames <- T.colnames %>% unique()
T.colnames <- T.colnames[-1]
total <- T.colnames %>% length()
Normalized.callus.data <- data.frame(AGI = callus.data$X,
                                     stringsAsFactors = F)
j <- 1
for (j in j:total) {
  Normalized.callus.data <- cbind(Normalized.callus.data, 
                                  callus.data %>% 
                                    select(starts_with(T.colnames[j])) %>% 
                                    apply(MARGIN = 1, FUN = mean) %>% 
                                    scale()
                                  )
  j <- j+1
}
colnames(Normalized.callus.data)[2:ncol(Normalized.callus.data)] <- T.colnames
#processing data----
T.data <- Normalized.Noshioka.data %>% 
  left_join(y = Normalized.callus.data, by = "AGI")
#PCA mean ver----
PCA.data <- T.data %>% select(-AGI)
PCA.data <- PCA.data %>% na.omit()
PCA.data <- PCA.data %>% lapply(FUN = as.numeric) %>% data.frame()
PCA.data <- PCA.data %>% t()
PCA <- prcomp(PCA.data, scale = TRUE)
T.component <- summary(PCA)
T.component <- T.component$importance[2, 1:2]
#ggplot----
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
ggsave(filename = "~/nishioka/Image/Nishioka&Callus_PCA.png", plot = g, width = 9, height = 7)
#heatmap----
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
ggsave(filename = "~/nishioka/Image/PCA_mean_heatmap.png", plot = g, width = 12, height = 18)