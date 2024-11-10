##conda activate celescope
setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/annotation/allele/TFBS")

library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(stringr)
snps.mb <- snps.from.file(file = "sig05_case_AD_SC_vid.bed",search.genome = BSgenome.Hsapiens.UCSC.hg19,format = "bed",indels = TRUE)

library(MotifDb)

#HOCOMOCOv11#
hocomo=subset(MotifDb, organism=='Hsapiens' & (dataSource == 'HOCOMOCOv11-core-A'|dataSource == 'HOCOMOCOv11-core-B'|dataSource == 'HOCOMOCOv11-core-C'|dataSource == 'HOCOMOCOv11-core-D'|dataSource == 'HOCOMOCOv11-secondary-A'|dataSource == 'HOCOMOCOv11-secondary-B'|dataSource == 'HOCOMOCOv11-secondary-C'|dataSource == 'HOCOMOCOv11-secondary-D'))

results_hocomo <- motifbreakR(snpList = snps.mb, filterp = TRUE,legacy = FALSE,
                       pwmList = hocomo,
                       threshold = 1e-4,
                       method = "default",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::MulticoreParam())
write.table(results_hocomo,"results_hocomo.txt",row.names=T,col.names=T,sep="\t",quote=T)

#plot#
pdf("chr4:89885714:T:C.pdf")
plotMB(results_hocomo, "chr4:89885714:T:C", reverseMotif = TRUE, effect = "strong")
dev.off()

results_hocomo$vid=paste(substring(results_hocomo$seqnames,4,nchar(results_hocomo$seqnames)),results_hocomo$start,results_hocomo$REF,results_hocomo$ALT,sep=":")

results_hocomo=unique(subset(results_hocomo,effect=="strong" & (dataSource == 'HOCOMOCOv11-core-A'|dataSource == 'HOCOMOCOv11-core-B'|dataSource == 'HOCOMOCOv11-core-C'|dataSource == 'HOCOMOCOv11-secondary-A'|dataSource == 'HOCOMOCOv11-secondary-B'|dataSource == 'HOCOMOCOv11-secondary-C'))[,c("vid","geneSymbol")])

TF=unlist(sapply(1:length(results_hocomo$vid),function(x){str_c(subset(results_hocomo,vid==results_hocomo$vid[x])$geneSymbol,collapse=";")}))
vid=results_hocomo$vid

TFBS_HOCOMOCO=unique(data.frame("vid"=vid,"TFBS_HOCOMOCO"=TF))
TFBS_HOCOMOCO=as.data.frame(TFBS_HOCOMOCO)

#jaspar2018#
jaspar=subset(MotifDb, organism=='Hsapiens' & dataSource == 'jaspar2018')

results_jaspar <- motifbreakR(snpList = snps.mb, filterp = TRUE,legacy = FALSE,
                       pwmList = jaspar,
                       threshold = 1e-4,
                       method = "default",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())
write.table(results_jaspar,"results_jaspar.txt",row.names=F,col.names=T,sep="\t",quote=F)

#plot#
pdf("4:89885714:T:C.pdf")
plotMB(results_jaspar, 4:89885714:T:C, reverseMotif = TRUE, effect = c("strong"))
dev.off()


results_jaspar$vid=paste(substring(results_jaspar$seqnames,4,nchar(results_jaspar$seqnames)),results_jaspar$start,results_jaspar$REF,results_jaspar$ALT,sep=":")

results_jaspar=unique(subset(results_jaspar,effect=="strong" )[,c("vid","geneSymbol")])

TF=unlist(sapply(1:length(results_jaspar$vid),function(x){str_c(subset(results_jaspar,vid==results_jaspar$vid[x])$geneSymbol,collapse=";")}))
vid=results_jaspar$vid

TFBS_jASPAR=unique(data.frame("vid"=vid,"TFBS_jASPAR"=TF))
TFBS_jASPAR=as.data.frame(TFBS_jASPAR)

TF=merge(TFBS_HOCOMOCO,TFBS_jASPAR,all=TRUE)
write.table(TF,"TFBS_motifbreakR.txt",row.names=F,col.names=T,sep="\t",quote=F)

save.image("TFBS_motifbreakR.Rdata")
