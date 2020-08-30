"~/nishioka/script/GO_analysis.R"
#load package----
library(dplyr)
library(stringr)
library(ggplot2)
#input data----
GO.Table <- readRDS(file = "~/bigdata/yasue/CY151620Exp_CoEXPNet/RDS/GOTable.rds")
#processing data----
GO.Table <- GO.Table %>% filter(M >= 10, FDR < 0.05)
T.Aspect <- GO.Table$Aspect %>% unique()
T.MCLNum <- GO.Table$MCLNum %>% unique()
i <- 1
for(i in i:length(T.MCLNum)){
  df <- c()
  j <- 1
  for(j in j:length(T.Aspect)){
    T.data <- GO.Table %>% filter(MCLNum == T.MCLNum[i], Aspect == T.Aspect[j]) %>% arrange(FDR) %>% select(Aspect, GOID, GOTerm, FDR)
    df <- rbind(df, T.data[1:5, ])
    j <- j+1
  }
  df <- df %>% na.omit()
  df <- df %>% mutate(EnrichmentScore = -log2(df$FDR))
  df$GOTerm <- factor(x = df$GOTerm, levels = df$GOTerm[order(df$EnrichmentScore, decreasing = T)], ordered = TRUE)
  g <- ggplot(df, aes(y = EnrichmentScore, x = GOTerm))
  g <- g + theme_bw()
  g <- g + geom_bar(stat = "identity")
  g <- g + coord_flip()
  g <- g + xlab("GO Term")
  g <- g + ylab("-log2(FDR)")
  g <- g + facet_grid(Aspect~., scales="free_y", space="free", labeller = label_both)
  g <- g + ggtitle(paste0("MCLNum", T.MCLNum[i]))
  #g <- g + scale_y_continuous(limits = c(0, 50))
  plot(g)
  title <- paste0("~/bigdata/yasue/CY151620Exp_CoEXPNet/Image/GO_Results/", "CY151620_MCLNum", T.MCLNum[i], "_BarGraph.png")
  ggsave(filename = title, plot = g)
  i <- i+1
}