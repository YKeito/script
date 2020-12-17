#load library----
library(tidyverse)
#Input published data & formatting data----
DEGs.Ave <- read.table("~/IMO/base/20201217_Aveppm_DEGs.txt", sep = "\t", header = T, quote = "")
callus.data <- read.table(file = "~/IMO/PublishedData/RNA-Seq_callus.txt", sep = "\t", header = T, stringsAsFactors = F)
callus.data <- callus.data %>% filter(X %in% DEGs.Ave$AGI)
#Convert expression data to Zscore
T.colnames <- DEGs.Ave %>% colnames()
T.colnames <- T.colnames[c(-1, -12)]
total <- T.colnames %>% length()
Normalized.Noshioka.data <- data.frame(AGI = DEGs.Ave$AGI,
                                       stringsAsFactors = F)
j <- 1
for (j in j:total) {
  Normalized.Noshioka.data <- cbind(Normalized.Noshioka.data, 
                                    DEGs.Ave %>% 
                                      select(starts_with(T.colnames[j])) %>% 
                                      scale()
  )
  j <- j+1
}
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
colnames(T.data) <- c("AGI", "Col_0d", "Col_2d", "Col_4d", "Col_7d", "Col_10d",
                      "rpt5a_0d", "rpt5a_2d", "rpt5a_4d", "rpt5a_7d", "rpt5a_10d",
                      "Col_C0", "Col_C14")
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
                 Group = rownames(PCA$x),
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
ggsave(filename = "~/IMO/Image/PCA/20201217Nishioka&Callus_Normalized_Aveppm_PCA.png", plot = g, width = 9, height = 7)
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
ggsave(filename = "~/IMO/Image/Heatmap/20201217Nishioka&Callus_Normalized_Aveppm_heatmap.png", plot = g, width = 12, height = 18)