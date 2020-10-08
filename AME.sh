
#FDR0.01
##iterations2.5
ls /home/yasue/AME_upstream500/DEG_FDR0.01/Motif/MultiFasta/Target/ > /home/yasue/AME_upstream500/DEG_FDR0.01/Motif/MultiFasta/allfasta.txt
for i in `cat /home/yasue/AME_upstream500/DEG_FDR0.01/Motif/MultiFasta/allfasta.txt`
do
ame --o /home/yasue/AME_upstream500/DEG_FDR0.01/Motif/Results/${i}_results --control /home/yasue/AME_upstream500/DEG_FDR0.01/Motif/MultiFasta/Control/control_${i} /home/yasue/AME_upstream500/DEG_FDR0.01/Motif/MultiFasta/Target/${i} /home/yasue/motif_databases/CIS-BP/Arabidopsis_thaliana.meme
done

#FDR0.05
##iterations2.5
ls /home/yasue/AME_upstream500/DEG_FDR0.05/inflations2.5/Motif/MultiFasta/Target/ > /home/yasue/AME_upstream500/DEG_FDR0.05/inflations2.5/Motif/MultiFasta/allfasta.txt
for i in `cat /home/yasue/AME_upstream500/DEG_FDR0.05/inflations2.5/Motif/MultiFasta/allfasta.txt`
do
ame --o /home/yasue/AME_upstream500/DEG_FDR0.05/inflations2.5/Motif/Results/${i}_results --control /home/yasue/AME_upstream500/DEG_FDR0.05/inflations2.5/Motif/MultiFasta/Control/control_${i} /home/yasue/AME_upstream500/DEG_FDR0.05/inflations2.5/Motif/MultiFasta/Target/${i} /home/yasue/motif_databases/CIS-BP/Arabidopsis_thaliana.meme
done
##iterations3.5
ls /home/yasue/AME_upstream500/DEG_FDR0.05/inflations3.5/Motif/MultiFasta/Target/ > /home/yasue/AME_upstream500/DEG_FDR0.05/inflations3.5/Motif/MultiFasta/allfasta.txt
for i in `cat /home/yasue/AME_upstream500/DEG_FDR0.05/inflations3.5/Motif/MultiFasta/allfasta.txt`
do
ame --o /home/yasue/AME_upstream500/DEG_FDR0.05/inflations3.5/Motif/Results/${i}_results --control /home/yasue/AME_upstream500/DEG_FDR0.05/inflations3.5/Motif/MultiFasta/Control/control_${i} /home/yasue/AME_upstream500/DEG_FDR0.05/inflations3.5/Motif/MultiFasta/Target/${i} /home/yasue/motif_databases/CIS-BP/Arabidopsis_thaliana.meme
done
##iterations4.5
ls /home/yasue/AME_upstream500/DEG_FDR0.05/inflations4.5/Motif/MultiFasta/Target/ > /home/yasue/AME_upstream500/DEG_FDR0.05/inflations4.5/Motif/MultiFasta/allfasta.txt
for i in `cat /home/yasue/AME_upstream500/DEG_FDR0.05/inflations4.5/Motif/MultiFasta/allfasta.txt`
do
ame --o /home/yasue/AME_upstream500/DEG_FDR0.05/inflations4.5/Motif/Results/${i}_results --control /home/yasue/AME_upstream500/DEG_FDR0.05/inflations4.5/Motif/MultiFasta/Control/control_${i} /home/yasue/AME_upstream500/DEG_FDR0.05/inflations4.5/Motif/MultiFasta/Target/${i} /home/yasue/motif_databases/CIS-BP/Arabidopsis_thaliana.meme
done
