####"C:/Users/Matsnaga-21/Desktop/R script (nishioka)/"####


#0.05
allRNASeq_cytoscape_th_0.05 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.05, ]
allRNASeq_cytoscape_th_possitive_0.05 <- allRNASeq_cytoscape_th_0.05[allRNASeq_cytoscape_th_0.05$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.05 <- allRNASeq_cytoscape_th_possitive_0.05[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.05, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.05.txt", append=F, quote = F, sep = "\t", row.names = F)

#0.01
allRNASeq_cytoscape_th_0.01 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.01, ]
allRNASeq_cytoscape_th_possitive_0.01 <- allRNASeq_cytoscape_th_0.01[allRNASeq_cytoscape_th_0.01$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.01 <- allRNASeq_cytoscape_th_possitive_0.01[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.01, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.01.txt", append=F, quote = F, sep = "\t", row.names = F)

#0.005
allRNASeq_cytoscape_th_0.005 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.005, ]
allRNASeq_cytoscape_th_possitive_0.005 <- allRNASeq_cytoscape_th_0.005[allRNASeq_cytoscape_th_0.005$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.005 <- allRNASeq_cytoscape_th_possitive_0.005[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.005, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.005.txt", append=F, quote = F, sep = "\t", row.names = F)

#0.001
allRNASeq_cytoscape_th_0.001 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.001, ]
allRNASeq_cytoscape_th_possitive_0.001 <- allRNASeq_cytoscape_th_0.001[allRNASeq_cytoscape_th_0.001$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.001 <- allRNASeq_cytoscape_th_possitive_0.001[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.001, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.001.txt", append=F, quote = F, sep = "\t", row.names = F)

#0.0005
allRNASeq_cytoscape_th_0.0005 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.0005, ]
allRNASeq_cytoscape_th_possitive_0.0005 <- allRNASeq_cytoscape_th_0.0005[allRNASeq_cytoscape_th_0.0005$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.0005 <- allRNASeq_cytoscape_th_possitive_0.0005[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.0005, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.0005.txt", append=F, quote = F, sep = "\t", row.names = F)

#0.0001
allRNASeq_cytoscape_th_0.0001 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.0001, ]
allRNASeq_cytoscape_th_possitive_0.0001 <- allRNASeq_cytoscape_th_0.0001[allRNASeq_cytoscape_th_0.0001$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.0001 <- allRNASeq_cytoscape_th_possitive_0.0001[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.0001, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.0001.txt", append=F, quote = F, sep = "\t", row.names = F)


#0.00005
allRNASeq_cytoscape_th_0.00005 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.00005, ]
allRNASeq_cytoscape_th_possitive_0.00005 <- allRNASeq_cytoscape_th_0.00005[allRNASeq_cytoscape_th_0.00005$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.00005 <- allRNASeq_cytoscape_th_possitive_0.00005[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.00005, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.00005.txt", append=F, quote = F, sep = "\t", row.names = F)

#0.00001
allRNASeq_cytoscape_th_0.00001 <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.00001, ]
allRNASeq_cytoscape_th_possitive_0.00001 <- allRNASeq_cytoscape_th_0.00001[allRNASeq_cytoscape_th_0.00001$interaction_value > 0, ]
exp_q_allRNASeq_cytoscape_possitive_th_0.00001 <- allRNASeq_cytoscape_th_possitive_0.00001[,-4:-5]
write.table(exp_q_allRNASeq_cytoscape_possitive_th_0.00001, file = "GRN_output/03_PCC calculation/200827coexptable_possitive_0.00001.txt", append=F, quote = F, sep = "\t", row.names = F)
