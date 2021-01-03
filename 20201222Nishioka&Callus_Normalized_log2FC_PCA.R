#load library----
library(tidyverse)
#Input published data & formatting data----
DEGs.Ave <- read.table("~/IMO/base/20201217_Aveppm_DEGs.txt", sep = "\t", header = T, quote = "")
Nishioka.data <- read.table("~/IMO/base/20201217_log2FC.txt", sep = "\t", header = T, quote = "")
Nishioka.data <- Nishioka.data %>% filter(reference_id %in% DEGs.Ave$AGI)
Nishioka.data <- Nishioka.data %>% select(-6:-18)
Normalized.Nishioka.data <- data.frame(AGI = DEGs.Ave$AGI,
                                       stringsAsFactors = F)
total <- ncol(Nishioka.data)
j <- 2
for (j in j:total) {
  Normalized.Nishioka.data <- cbind(Normalized.Nishioka.data, 
                                    Nishioka.data %>% 
                                      select(j) %>% 
                                      scale()
  )
  j <- j+1
}
colnames(Normalized.Nishioka.data)[2:ncol(Nishioka.data)] <- colnames(Nishioka.data)[2:ncol(Nishioka.data)]
#callus-----
callus.data <- read.table(file = "~/IMO/PublishedData/RNA-Seq_callus.txt", sep = "\t", header = T, stringsAsFactors = F)
callus.data <- callus.data %>% filter(X %in% DEGs.Ave$AGI)

C0.data <- callus.data %>% 
  select(starts_with("Col_C0")) %>% 
  apply(MARGIN = 1, FUN = mean)
C0.data <- C0.data + 0.001
C14.data <- callus.data %>% 
  select(starts_with("Col_C14")) %>% 
  apply(MARGIN = 1, FUN = mean)

C14.data <- C14.data + 0.001
ratio <- C14.data/C0.data

Normalized.callus.data <- data.frame(AGI = callus.data$X,
                                     X0_14d_log2FC = log2(ratio),
                                     stringsAsFactors = F)
Normalized.callus.data$X0_14d_log2FC <- Normalized.callus.data$X0_14d_log2FC %>% scale()
#join--------
T.data <- Normalized.Nishioka.data %>% 
  left_join(y = Normalized.callus.data, by = "AGI")
#PCA----
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
df <- df %>% mutate(Sample = c("rpt5a", "rpt5a", "rpt5a", "rpt5a", "callus"))
df <- df %>% mutate(Time = str_split(df$Group, pattern = "_log2FC", simplify = T)[, 1])
g <- ggplot(data = df, aes(x=PC1, y=PC2, color = Time, label = Group, shape = Sample))
g <- g + theme_bw()
g <- g + geom_point(size = 5)
g <- g + xlab(paste0("PC1 (", formatC(T.component[1]*100, width = 3), "%)"))
g <- g + ylab(paste0("PC2 (", formatC(T.component[2]*100, width = 3), "%)"))
g <- g + theme(axis.text = element_text(size=18))
g <- g + theme(axis.title = element_text(size=18))
ggsave(filename = "~/IMO/Image/PCA/20201222Nishioka&Callus_Normalized_log2FC_PCA.png", plot = g, width = 9, height = 7)
#heatmap----
colnames(T.data)[2:ncol(T.data)] <- paste0(c(rep("rpt5a_", 4), "callus_"), colnames(T.data)[2:ncol(T.data)])
colnames(T.data)[2:ncol(T.data)] <- str_split(colnames(T.data)[2:ncol(T.data)], pattern = "_log2FC", simplify = T)[, 1]
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
ggsave(filename = "~/IMO/Image/Heatmap/20201222Nishioka&Callus_Normalized_log2FC_heatmap.png", plot = g, width = 12, height = 18)