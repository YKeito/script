library(tidyverse)
T.FDR <- c("FDR005", "FDR001")
path <- c("~/IMO/RDS/200909(FDR0.05)PCC_CytoscapeFormat.rds", "~/IMO/RDS/200909(FDR0.01)PCC_CytoscapeFormat.rds")
i <- 1
for (i in i:length(path)) {
  T.data <- readRDS(file = path[i])
  T.data <- T.data %>% filter(interaction_value > 0)
  th.FDR <- c(5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5)
  Numedges <- c()
  Numnodes <- c()
  j <- 1
  for (j in j:length(th.FDR)) {
    th.data <- T.data %>% filter(q_value < th.FDR[j])
    Numedges <- c(Numedges, nrow(th.data))
    Numnodes <- c(Numnodes, length(union(T.data$source_genes, T.data$target_genes)))
    j <- j + 1
  }
  ggplot.data <- data.frame(th = th.FDR,
                            NumEdges = Numedges,
                            stringsAsFactors = F
  )
  g <- ggplot(ggplot.data, 
              aes(x = th, y = NumEdges)
  )
  g <- g + geom_point(size = 3)
  g <- g + theme_bw()
  g <- g + geom_hline(yintercept = 6e5, 
                      colour="red", 
                      size=1.5, 
                      linetype=2
  )
  g <- g + labs(title = "PCC > 0 and Number of edges in various thresholds of FDR",
                x = "PCC > 0 and FDR used for cut-off",
                y = "Numbers of edges"
  )
  g <- g + scale_y_continuous(breaks=seq(0,max(Numedges),5e5))
  plot(g)
  title <- paste0("~/IMO/Image/DeterminationOfThreshold/", T.FDR[i], "_NumEdges_FDR_CutOff.png")
  ggsave(filename = title)
  print(Numnodes)
  i <- i + 1
}
T.data <- readRDS(file = "~/IMO/RDS/200909(FDR0.05)PCC_CytoscapeFormat.rds")
T.data <- T.data %>% filter(interaction_value > 0, q_value < 1e-04)
write.table(x = T.data,
            file = "~/IMO/Table/PCC_FDR1e-4&PCC_Possitive.txt", sep = "\t", quote = F, row.names = F)
T.data <- readRDS(file = "~/IMO/RDS/200909(FDR0.01)PCC_CytoscapeFormat.rds")
T.data <- T.data %>% filter(interaction_value > 0, q_value < 1e-03)
write.table(x = T.data,
            file = "~/IMO/Table/PCC_FDR1e-3&PCC_Possitive.txt", sep = "\t", quote = F, row.names = F)
