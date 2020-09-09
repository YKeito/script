####pvalue&qvalue calculate####
####"C:/Users/Matsnaga-21/Desktop/R script (nishioka)/"####


####read.table####
ppm_Colvsrpt<- read.table("GRN_input/01_pqvalue_vsday0/200817iMO2d_vs0d.txt", sep = "\t", row.names = 1, header = T,quote = "")
x<-ppm_Colvsrpt[,1:3]
y<-ppm_Colvsrpt[,4:6]

#####pvalue calculate####
n <- 1
total <- nrow(ppm_Colvsrpt)
p_ans <- c()      ###c()�̓x�N�g��
for(n in n:total){
  out1<-t.test(x[n, ],y[n, ], var.equal=T)
  #out1$p.value
  #print(out1$p.value)
  p_ans[n]<-out1$p.value
  n <- n+1
}

##### NA are replaced by 1 ####
p_ans[is.na(p_ans)]<-1

#####qvalue calculate (p.adjust) ####
q_ans<-p.adjust(p_ans,method="BH")


####p value��q value �f�[�^�̌���#####
gene <- ppm_Colvsrpt[,0]
pq_ans <- data.frame(gene,p_ans,q_ans)

####output####
write.table(pq_ans,"GRN_output/01_pqvalue_vsday0/190214day10pqvalue.txt",quote=F,row.names=T,sep="\t", col.names=T)
