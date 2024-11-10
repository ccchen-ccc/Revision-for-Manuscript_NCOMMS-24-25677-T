#####fig3h
#########copd data from hum0014.v17.COPD.v1 https://humandbs.biosciencedbc.jp/files/hum0014/hum0014.v17.COPD.v1.zip
#awk 'NR==1 || ($1==4 && $2>88000001 && $2<93700000)' COPD.auto.rsq07.mac10.txt >COPD.bbj_4q22.txt ###
##########
setwd("/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step3/data")
library(tidyverse)
library(ggplot2)
lc <- fread('/public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3h/asso_combined_6_rpWGS_chr_all_case.firth.covar1_source7_GC1.1.adjusted.filter.output')
lc <- lc %>% 
  separate(ID,into = c('chr','pos','ref','alt'),remove = F) %>% 
  filter(chr==4,pos>=88000001&pos<=93700000)
write.table(lc,"4q22.1.asso",row.names=F,col.names=T,quote=F,sep="\t")

bbj <- read.table('COPD.bbj_4q22.txt',header = T,sep = ' ')[c(1,2,4,5,9,10,12)]
colnames(bbj)[c(5:7)] <-paste0('bbj_copd_',colnames(bbj)[c(5:7)])
bbj$CHR=as.character(bbj$CHR)
bbj$POS=as.character(bbj$POS)
test <- lc %>% 
  left_join(bbj,by = c('chr'='CHR','pos'='POS','ref'='Allele1','alt'='Allele2')) %>% 
  filter(GC<=0.05)
cor_p5 = ggplot(data = test, aes(x = -log10(GC), y = -log10(bbj_copd_p.value))) + 
  geom_point() + geom_smooth(method = glm) + theme_classic2()+
  scale_color_manual(values = c('#FF7400')) + 
  xlab("LC_GWAS_log10(p)") + ylab("BBJ_COPD_log10(p)")+
  stat_cor(method = 'pearson', aes(x = -log10(GC), y = -log10(bbj_copd_p.value))) 
ggsave(cor_p5,filename = 'fig3h.pdf',width = 10,height = 10,dpi = 300,path = './')