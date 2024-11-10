setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data/LungENN/prediction")
library(data.table)
dat1=fread("sig05_case_AD_SC_diffs_v1.tsv")
dat1=as.data.frame(dat1)
dat1$max=sapply(1:nrow(dat1),function(x){max(dat1[x,9:298])})
dat1$min=sapply(1:nrow(dat1),function(x){min(dat1[x,9:298])})

dat1$LungENN_v1_all_max=ifelse(abs(dat1$max)>abs(dat1$min),dat1$max,dat1$min)

dat1=dat1[,c("name","LungENN_v1_all_max")]

write.table(dat1,"../sig05_case_AD_SC_vid_LungENN.txt",row.names=F,col.names=T,sep="\t",quote=F)
