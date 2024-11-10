library(tidyverse)
library(data.table)
library(patchwork)
library(ggcorrplot)
library(ggpubr)
library(vroom)
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

oligoCounts = fread(paste(args[2],"_Oligo_counts_sum.txt",sep=""),header=T)
oligoCounts = oligoCounts[which(rowMeans(oligoCounts[,c(5:16,29:40,53:64)]) > 150),]

mergedOligoCPM = oligoCounts
mergedOligoCPM[is.na(mergedOligoCPM)] <- 0

mergedOligoCPM[,c(5:16,29:40,53:64)] =sweep(mergedOligoCPM[,c(5:16,29:40,53:64)],2,colSums(mergedOligoCPM[,c(5:16,29:40,53:64)])/1000000,`/`)

dna1vdna2=ggplot(data = mergedOligoCPM, aes(x = DNA1_counts_ref, y = DNA2_counts_ref)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("DNA Replicate 1 (Reference Allele)") + ylab("DNA Replicate 2 (Reference Allele)")+
				scale_x_continuous(limit = c(0,1000)) + scale_y_continuous(limit = c(0,1000)) +
				stat_cor(method = 'pearson', aes(x = DNA1_counts_ref, y = DNA2_counts_ref))
rna1vrna2=ggplot(data = mergedOligoCPM, aes(x = cDNA1_counts_ref, y = cDNA2_counts_ref)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("cDNA Replicate 1 (Reference Allele)") + ylab("cDNA Replicate 2 (Reference Allele)")+
				scale_x_continuous(limit = c(0,1000)) + scale_y_continuous(limit = c(0,1000)) +
				stat_cor(method = 'pearson', aes(x = cDNA1_counts_ref, y = cDNA2_counts_ref))
dna1vrna1_ref=ggplot(data = mergedOligoCPM, aes(x = DNA1_counts_ref, y = cDNA1_counts_ref)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("DNA Replicate 1 (Reference Allele)") + ylab("cDNA Replicate 1 (Reference Allele)")+
				scale_x_continuous(limit = c(0,1000)) + scale_y_continuous(limit = c(0,1000)) +
				stat_cor(method = 'pearson', aes(x = DNA1_counts_ref, y = cDNA1_counts_ref))
dna1vrna1_alt=ggplot(data = mergedOligoCPM, aes(x = DNA1_counts_alt, y = cDNA1_counts_alt)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("DNA Replicate 1 (Reference Allele)") + ylab("cDNA Replicate 1 (Reference Allele)")+
				scale_x_continuous(limit = c(0,1000)) + scale_y_continuous(limit = c(0,1000)) +
				stat_cor(method = 'pearson', aes(x = DNA1_counts_alt, y = cDNA1_counts_alt))
dnaRefAlt=ggplot(data = mergedOligoCPM, aes(x = DNA1_counts_ref, y = DNA1_counts_alt)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("DNA Replicate 1 (Reference Allele)") + ylab("DNA Replicate 1 (Alternative Allele)")+
				scale_x_continuous(limit = c(0,1000)) + scale_y_continuous(limit = c(0,1000)) +
				stat_cor(method = 'pearson', aes(x = DNA1_counts_ref, y = DNA1_counts_alt))
rnaRefAlt=ggplot(data = mergedOligoCPM, aes(x = cDNA1_counts_ref, y = cDNA1_counts_alt)) + 
				geom_point() + geom_smooth(method = lm) + theme_bw()+
				scale_color_manual(values = c('#FF7400')) + 
				xlab("cDNA Replicate 1 (Reference Allele)") + ylab("cDNA Replicate 1 (Alternative Allele)")+
				scale_x_continuous(limit = c(0,1000)) + scale_y_continuous(limit = c(0,1000)) +
				stat_cor(method = 'pearson', aes(x = cDNA1_counts_ref, y = cDNA1_counts_alt))

ggsave(paste(args[2],"_dna1vdna2.bmp",sep=""),dna1vdna2,dpi=300)
ggsave(paste(args[2],"_dna1vrna1_ref.bmp",sep=""),dna1vrna1_ref,dpi=300)
ggsave(paste(args[2],"_dna1vrna1_alt.bmp",sep=""),dna1vrna1_alt,dpi=300)
ggsave(paste(args[2],"_dnaRefAlt.bmp",sep=""),dnaRefAlt,dpi=300)
ggsave(paste(args[2],"_rnaRefAlt.bmp",sep=""),rnaRefAlt,dpi=300)
ggsave(paste(args[2],"_rna1vrna2.bmp",sep=""),rna1vrna2,dpi=300)

corMat = cor(mergedOligoCPM[,c(5:16,29:40)], method = "pearson", use = "pairwise.complete.obs")
samples = c("DNA 1 Ref", "DNA 2 Ref","DNA 3 Ref","DNA 4 Ref", "DNA 5 Ref","DNA 6 Ref",
			"cDNA 1 Ref","cDNA 2 Ref","cDNA 3 Ref","cDNA 4 Ref","cDNA 5 Ref","cDNA 6 Ref",
			"DNA 1 Alt", "DNA 2 Alt","DNA 3 Alt","DNA 4 Alt", "DNA 5 Alt","DNA 6 Alt",
			"cDNA 1 Alt","cDNA 2 Alt","cDNA 3 Alt","cDNA 4 Alt","cDNA 5 Alt","cDNA 6 Alt")

colnames(corMat) = samples
rownames(corMat) = samples

options(repr.plot.width = 10, repr.plot.height = 10)
corrplot=ggcorrplot(corMat,
           hc.order = F,
           tl.cex = 10,
           pch.cex = 3,
           lab_size = 2,
           lab = T,
           show.diag = T,
           type = "upper") + 
    theme_pubr(x.text.angle = 45, base_size = 10) +
    scale_fill_gradient2(limit = c(0,1), midpoint = 0.5,low = "blue",mid = "white", high =  "red", name = "Pearson Correlation   Coefficient") +
    theme(legend.position = c(0.75, 0.3)) +
    xlab("") + ylab("")
ggsave(paste(args[2],"_corrplot.bmp",sep=""),corrplot,dpi=300)
