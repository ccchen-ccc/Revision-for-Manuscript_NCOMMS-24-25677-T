library(data.table)
library(reshape)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

setwd(paste("/public/home/ccc2018/data/project/MPRA/LC/GWAS/annotation/allele/fvars_LC/eqtl/boxplot/",args[1],"/plink",sep=""))
dat1=fread("/public/home/ccc2018/data/database/LC_WGS_germline/eqtl_v2/normal/lung_normal_expression.bed.gz")
dat1=as.data.frame(dat1)
dat1=dat1[,-c(1:3)]
row.names(dat1)=dat1[,1]
dat2=t(dat1[,-1])
dat2<-transform(dat2,FID = rownames(dat2))


pair=fread(paste0("./",args[2],".lung_normal_sigeqtl.txt",sep=""))
anno=fread("/public/home/ccc2018/data/database/ref/gencodev19/genecodev19_anno.txt")
anno=subset(anno,gene_id%in%pair$V2)[,c("symbol","gene_id")]

pair=merge(pair,anno,by.x="V2",by.y="gene_id")

for (i in 1:nrow(pair)){
newname<-c(pair$V2[i],"FID")
print (newname)
dat2_1<-dat2[newname]
dat3=fread("/public/home/ccc2018/data/database/LC_WGS_germline/eqtl_v2/normal/covar_normal.txt.gz")
dat3=as.data.frame(dat3)
row.names(dat3)=dat3[,1]
dat4=t(dat3[,-1])
dat4<-transform(dat4,FID = rownames(dat4))
rownames(dat4)=NULL
dat5=merge(dat2_1,dat4)
dat6=fread(paste0("./",args[2],".raw",sep=""))
#dat6$FID<-paste(sapply(strsplit(dat6$FID, '-'), '[', 1), sapply(strsplit(dat6$FID, '-'), '[' ,2), sep ='-')
dat7=merge(dat6,dat5)
#dat7=rename(dat7,c(get(pair$V2[i])=pair$symbol[i]))
dat7=as.data.frame(dat7)
colnames(dat7)

dat8=dat7[complete.cases(dat7[,7]),]
#r1<-lm(dat8[,8]~dat8[,9]+PC1+PC2+PC3+PC4+PC5+InferredCov1+InferredCov2+InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8+InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14+InferredCov15+gender,data=dat8)
#summary(r1)

plot=ggplot(dat8, aes(x = factor(dat8[,7]), y = dat8[,8], color=factor(dat8[,7]))) + 
  scale_x_discrete(labels=c("0","1",paste(sapply(strsplit(names(dat8)[7],"_"),'[',2),sapply(strsplit(names(dat8)[7],"_"),'[',2),sep=""))) + 
  geom_jitter(position=position_jitter(width=.2, height=0))+geom_boxplot(fill=NA, size=.9)+theme_bw()+
  annotate("text", family="Helvetica",x=2.6, y=3.4, label=paste("Beta=",round(pair$V8[i],2)," ",sep=""), hjust=0,size=3.5) +
  annotate("text", family="Helvetica",x=2.6, y=3.1, label=paste("P=",pair$V7[i]," ",sep=""), hjust=0,size=3.5) +
  theme(
    panel.background =  element_rect(fill = NA),
    legend.position="none",axis.text.y=element_text(size=12,color="#030303",face="bold"),
    axis.title.y=element_text(size=10,face="bold"),axis.text.x=element_text(size=12,color="#030303",face="bold"),
    panel.grid.major=element_line(colour=NA),panel.grid=element_blank())+
	ylab(paste("Normalized Expression of ",pair$symbol[i],sep=""))+xlab(NULL)+scale_colour_manual(values=c("#006699","#458B74","#DDA520"))
ggsave(file=paste(pair$symbol[i],"_",args[2],".pdf",sep=""),plot,dpi=300,width=4,height=3)
}
