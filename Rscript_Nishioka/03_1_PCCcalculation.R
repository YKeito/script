#PCCcalculation#
#02_DEGs (FDRcut).Rの続きで計算(PCC input data need only p value)#

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


####PCC qvalue cut####

allRNASeq_cytoscape_th <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.05, ]
allRNASeq_cytoscape_th_possitive <- allRNASeq_cytoscape_th[allRNASeq_cytoscape_th$interaction_value > 0, ]
#write.table(allRNASeq_cytoscape_th, file = "GRN_output/03_PCC calculation/200827_0.05_coexptable.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(allRNASeq_cytoscape_th_possitive, file = "GRN_output/03_PCC calculation/200827_0.05_coexptable_possitive.txt", append=F, quote = F, sep = "\t", row.names = F)