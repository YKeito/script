library(tidyverse)
library(gplots)

Nishioka.data <- read.table("~/IMO/base/20201217_log2FC.txt", sep = "\t", header = T, quote = "")
Nishioka.data <- Nishioka.data %>% select(-contains("pvalue"))
T.names <- c("2d", "4d", "7d", "10d")
T.label <- c("02d", "04d", "07d", "10d")
df <- c()
DEGs.list <- c()
i <- 1
for(i in i:length(T.names)){
  Nishioka.Exp <- Nishioka.data %>% select(contains(paste0(T.names[i], "_log2")))
  FDR.data <- Nishioka.data %>% select(contains(paste0(T.names[i], "_qvalue")))
  FDR.logic <- FDR.data < 0.05
  FDR.logic <- FDR.logic %>% apply(1, FUN = sum) 
  Nishioka.DEGs <- Nishioka.Exp %>% subset(FDR.logic == 2)
  
  Up.logic <- Nishioka.DEGs > 0
  Down.logic <- Nishioka.DEGs < 0
  T.DEGs.list <- Nishioka.data %>% select(reference_id) %>% subset(FDR.logic == 2)
  DEGs.list <- c(DEGs.list, list(T.DEGs.list$reference_id))
  df <- rbind(df,
              data.frame(Sample = T.label[i],
                         Group = c("02Down", "01Up"),
                         value = c(sum(Down.logic), sum(Up.logic)))
  )
  i <- i+1
}

g <- ggplot(df, aes(x = Sample, y = value, fill = Group))
g <- g + geom_bar(stat = "identity")
g <- g + theme_bw()
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) #軸名消す
g <- g + theme(legend.position = 'none') #凡例消す
g <- g + theme(axis.text = element_text(size = 18))
plot(g)
ggsave(filename = "~/IMO/Image/DEGs_Graph/20201217NumDEGs_M0VSMx&MxVSWx.png", plot = g, width = 8, height = 6)

#VennDiagram----
names(DEGs.list) <- T.label
Venn <- venn(DEGs.list)
names(attr(Venn,"intersections")) <- paste0(names(attr(Venn,"intersections")), ",")
T.Venn <- data.frame(GeneList = str_split(unlist(attr(Venn,"intersections"), use.names = F), pattern = ",", simplify = T),
                     conditions = str_split(names(unlist(attr(Venn,"intersections"))), pattern = ",", simplify = T)[, 1],
                     stringsAsFactors = F)
write.table(x = T.Venn, file = "~/IMO/Table/20201217Nishioka_VennDiagram.txt", sep = "\t", quote = F, row.names = F)
ggsave(filename = "~/IMO/Image/DEGs_Graph/20201217Nishioka_VennDiagram.png", plot = plot(Venn), width = 6, height = 4)
