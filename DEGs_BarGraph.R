#local R-studio version
#"~/nishioka/script/DEGs_BarGraph.R"
Nishioka.data <- read.table("~/nishioka/base/20200317_log2FC&FDR.txt", sep = "\t", header = T)
T.names <- c("0_2d", "0_4d", "0_7d", "0_10d")
df <- c()
i <- 1
for(i in i:length(T.names)){
  Nishioka.DEGs <- Nishioka.data %>% select(contains(T.names[i]))
  FDR.data <- Nishioka.DEGs %>% select(ends_with("qvalue"))
  FDR.logic <- FDR.data < 0.05
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
df$Sample <- factor(x = df$Sample,
                    levels = df$Sample[1:8],
                    ordered = TRUE)

g <- ggplot(df, aes(x = Sample, y = value, fill = Group))
g <- g + geom_bar(stat = "identity")
g <- g + theme_bw()
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) #Ž²–¼Á‚·
g <- g + theme(legend.position = 'none') #–}—áÁ‚·
g <- g + theme(axis.text = element_text(size = 18))
plot(g)
ggsave(filename = "~/nishioka/Image/NumDEGs_vs0_BarGraph.png", plot = g, width = 8, height = 6)