setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data")
library(data.table)
library(parallel)

anno=fread("sig05_case_AD_SC_60bp_expEffect_3cell_cov_rep6_0.01_filter.txt")
dat3=anno[order(-anno$padj_expr_tag),]

dat3$EncodeLung_peak_tag=ifelse(is.na(dat3$EncodeLung_peak)|dat3$EncodeLung_peak=="",0,1)

dat4=dat3[!duplicated(dat3$vid),c(1,11:67)]

func=function(z){
p_adj_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_tag==0&get(z)==0))),ncol=2))$p
OR_adj_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_tag==0&get(z)==0))),ncol=2))$estimate
ci95_low_adj_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_tag==0&get(z)==0))),ncol=2))$conf[1]
ci95_high_adj_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_tag==0&get(z)==0))),ncol=2))$conf[2]

p_adj_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==0))),ncol=2))$p
OR_adj_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==0))),ncol=2))$estimate
ci95_low_adj_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==0))),ncol=2))$conf[1]
ci95_high_adj_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_pos_tag==0&get(z)==0))),ncol=2))$conf[2]

p_adj_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==0))),ncol=2))$p
OR_adj_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==0))),ncol=2))$estimate
ci95_low_adj_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==0))),ncol=2))$conf[1]
ci95_high_adj_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,padj_expr_neg_tag==0&get(z)==0))),ncol=2))$conf[2]

p_value_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==0))),ncol=2))$p
OR_value_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==0))),ncol=2))$estimate
ci95_low_value_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==0))),ncol=2))$conf[1]
ci95_high_value_expr_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_tag==0&get(z)==0))),ncol=2))$conf[2]

p_value_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==0))),ncol=2))$p
OR_value_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==0))),ncol=2))$estimate
ci95_low_value_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==0))),ncol=2))$conf[1]
ci95_high_value_expr_pos_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_pos_tag==0&get(z)==0))),ncol=2))$conf[2]

p_value_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==0))),ncol=2))$p
OR_value_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==0))),ncol=2))$estimate
ci95_low_value_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==0))),ncol=2))$conf[1]
ci95_high_value_expr_neg_tag=fisher.test(matrix(c(nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==1&get(z)==0)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==1)),nrow(subset(dat4,pvalue_expr_neg_tag==0&get(z)==0))),ncol=2))$conf[2]

enrich=data.frame(tag=z,p_adj_expr_tag=p_adj_expr_tag,OR_adj_expr_tag=paste(round(OR_adj_expr_tag,2)," (",round(ci95_low_adj_expr_tag,2),"-",round(ci95_high_adj_expr_tag,2),")",sep=""),p_adj_expr_pos_tag=p_adj_expr_pos_tag,OR_adj_expr_pos_tag=paste(round(OR_adj_expr_pos_tag,2)," (",round(ci95_low_adj_expr_pos_tag,2),"-",round(ci95_high_adj_expr_pos_tag,2),")",sep=""),p_adj_expr_neg_tag=p_adj_expr_neg_tag,OR_adj_expr_neg_tag=paste(round(OR_adj_expr_neg_tag,2)," (",round(ci95_low_adj_expr_neg_tag,2),"-",round(ci95_high_adj_expr_neg_tag,2),")",sep=""),p_value_expr_tag=p_value_expr_tag,OR_value_expr_tag=paste(round(OR_value_expr_tag,2)," (",round(ci95_low_value_expr_tag,2),"-",round(ci95_high_value_expr_tag,2),")",sep=""),p_value_expr_pos_tag=p_value_expr_pos_tag,OR_value_expr_pos_tag=paste(round(OR_value_expr_pos_tag,2)," (",round(ci95_low_value_expr_pos_tag,2),"-",round(ci95_high_value_expr_pos_tag,2),")",sep=""),p_value_expr_neg_tag=p_value_expr_neg_tag,OR_value_expr_neg_tag=paste(round(OR_value_expr_neg_tag,2)," (",round(ci95_low_value_expr_neg_tag,2),"-",round(ci95_high_value_expr_neg_tag,2),")",sep=""))
return(enrich)
}
z <- names(dat4)[c(17:ncol(dat4))]
cl <- makeCluster(10,type="FORK")
res <- parLapply(cl,z,func)
results <- do.call('rbind',res)
stopCluster(cl)
write.table(results,"enrich_exp_3cell_cov_rep6_0.01_filter.txt",row.names=F,col.names=T,sep="\t",quote=F)
