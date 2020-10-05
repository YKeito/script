#load package----
library(tidyverse)
library(Hmisc)
#input base data----
allRNASeq <- subdata %>% select(1:4)
t <- proc.time()
CY.Exp <- t(CY.Exp)
base <- rcorr(as.matrix(CY.Exp), type = "pearson")
temp <- combn(colnames(CY.Exp), 2)
CY.PCC <- data.frame(source_genes = temp[1, ],
                     interaction_value = base$r[lower.tri(base$r)],
                     target_genes = temp[2, ],
                     p_value = base$P[lower.tri(base$P)],
                     q_value = p.adjust(base$P[lower.tri(base$P)], method = "BH"),
                     stringsAsFactors = F
)
saveRDS(object = CY.PCC, file = "~/bigdata/yasue/CY151620Exp_CoEXPNet/RDS//CY151620PCC.rds")
write.table(x = CY.PCC %>% select(-p_value) %>% filter(q_value < 5e-2),
            file = "~/bigdata/yasue/CY151620Exp_CoEXPNet/Table/CY151620PCC_FDR5e-2.txt", sep = "\t", quote = F, row.names = F)
#elapsed time----
after <- proc.time()
print(after - before)#175.43 sec

"""
#PCCcalculation#
#02_DEGs (FDRcut).Rの続きで計�?(PCC input data need only p value)#

####名称変更####
allRNASeq <- subdata

####進行状況確認用####
t <- proc.time()

####PCC calculation####
allRNASeq_foldchange <- allRNASeq[, 1:8]
test <- allRNASeq_foldchange[1:8, ]
n <- 1
m <- 2
base <- c()
PCC <- c()
PCC_all <- list()
PCC_pvalue <- c()
PCC_pvalue_all <- list()
total <- nrow(allRNASeq_foldchange)
for(n in n:c(total-1)){
  for(m in m:total){
    base <- cor.test(as.numeric(allRNASeq_foldchange[n, ]), as.numeric(allRNASeq_foldchange[m, ]), method = "pearson")
    PCC <- c(PCC, base$estimate)
    PCC_pvalue <- c(PCC_pvalue, base$p.value)
    m <- m+1
  }
  PCC_all <- c(PCC_all, list(PCC))
  PCC <- c()
  PCC_pvalue_all <- c(PCC_pvalue_all, list(PCC_pvalue))
  PCC_pvalue <- c()
  print(total-n)
  n <- n+1
  m <- n+1
}
t1 <- proc.time()-t
print(t1)
save.image(file = "GRN_output/03_PCC calculation/200827PCC(FDR0.05).RData")

####PCC qvalue calculation####
PCC_qvalue_all <- p.adjust(unlist(PCC_pvalue_all), method = "BH")

####Cytoscape####
test <- combn(rownames(allRNASeq_foldchange), 2)
source_genes <- test[1, ]
target_genes <- test[2, ]

allRNASeq_foldchange_cytoscape <- data.frame(source_genes = source_genes, 
                                             interaction_value = unlist(PCC_all), 
                                             target_genes = target_genes,
                                             p_value = unlist(PCC_pvalue_all),
                                             q_value = PCC_qvalue_all
)

####save####
saveRDS(allRNASeq_foldchange_cytoscape, file = "200909(FDR0.01)PCC_CytoscapeFormat.rds")

allRNASeq_foldchange_cytoscape0.05 <- readRDS("200909(FDR0.01)PCC_CytoscapeFormat.rds")

####PCC qvalue cut####
allRNASeq_cytoscape_th <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.05, ]
allRNASeq_cytoscape_th_possitive <- allRNASeq_cytoscape_th[allRNASeq_cytoscape_th$interaction_value > 0, ]
#write.table(allRNASeq_cytoscape_th, file = "GRN_output/03_PCC calculation/200827_0.05_coexptable.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(allRNASeq_cytoscape_th_possitive, file = "GRN_output/03_PCC calculation/200827_0.05_coexptable_possitive.txt", append=F, quote = F, sep = "\t", row.names = F)
"""