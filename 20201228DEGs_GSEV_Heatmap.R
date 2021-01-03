#Input data----
DEGs.Enrichment <- readRDS(file = "~/IMO/RDS/20201228DEGs_Enrichment.rds")
#processing data----
T.data <- DEGs.Enrichment %>% 
  filter(M >= 10, n >= 10) %>% 
  select(- Intersection_AGI) %>% 
  mutate(count = ifelse(pvalue <= 0.05, 1, 0))

T.DEGs.summary <- T.data %>% 
  group_by(day) %>% 
  summarise(Summary = sum(count)) %>% 
  filter(Summary != 0)

T.MCL.summary <- T.data %>% 
  group_by(MCLNum) %>% 
  summarise(Summary = sum(count)) %>% 
  filter(Summary != 0)

T.data <- T.data %>% 
  filter(day %in% T.DEGs.summary$day, MCLNum %in% T.MCL.summary$MCLNum)

df <- T.data %>% select(day, MCLNum, pvalue)
T.matrix <- df %>% spread(key = MCLNum, value = pvalue)
rownames(T.matrix) <- T.matrix$day
T.matrix <- T.matrix %>% select(-day)
res <- hclust(dist(T.matrix), method = "ward.D2")
df$day <- factor(x = df$day,
                        levels = res$labels[res$order],
                        ordered = TRUE)
res.MCL <- hclust(dist(t(T.matrix)), method = "ward.D2")
df$MCLNum <- factor(x = df$MCLNum,
                    levels = res.MCL$labels[res.MCL$order],
                    ordered = TRUE)
df$pvalue[df$pvalue > 5e-2] <- 1
df$pvalue[df$pvalue <= 5e-2 & df$pvalue > 5e-5] <- 0.67
df$pvalue[df$pvalue <= 5e-5 & df$pvalue > 5e-10] <- 0.34
df$pvalue[df$pvalue <= 5e-10] <- 0
T.MCLNum <- df$MCLNum %>% unique()
T.data <- c()
i <- 1
for(i in i:length(T.MCLNum)){
  test <- sum(df$pvalue[df$MCLNum == T.MCLNum[i]] == 1) != df$day %>% unique() %>% length()
  if(test){
    T.data <- rbind(T.data,
                    df %>% filter(MCLNum == T.MCLNum[i])
    )
  }
  i <- i+1
}
g <- ggplot(data = T.data, aes(x = MCLNum, y = day, fill = pvalue))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient(high="white",low="red", na.value = "white")
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g <- g + theme(axis.ticks.y = element_blank())
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + theme(legend.title = element_blank())
g <- g + theme(legend.text = element_text(size=7))
g <- g + coord_flip()
plot(g)
ggsave(filename = "~/IMO/Image/20201228DEGs_GSEV.png", plot = g, width = 1, height = 6)