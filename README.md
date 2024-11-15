# Massively-parallel-variant-to-function-mapping-determines-regulatory-functional-variants-of-NSCLC
Welcome to the GitHub repository for the project titled "Massively parallel variant-to-function mapping determines regulatory functional variants of non-small cell lung cancer." This project focuses on identifying and understanding the regulatory functional variants associated with non-small cell lung cancer (NSCLC) using MPRA and advanced computational techniques.

Data and code for files are organized in the /manuscript directory.

######Step 1 MPRA data analysis
####Analysis of sequencing data from MPRA experiments mainly drew on previous study [PMID- 35298243, https://www.science.org/doi/10.1126/science.abj5117], and the code was avalible from "10.5281/zenodo.5921041"
####MPRA results are presented in /manuscript/step1/data/vid_alleleEffect_for_MPRA.txt

######Step 2. frVars identification and annotations
####Lung-related annotations:DNase-seq, ATAC-seq, active histone ChIP-Seq and TF-ChIP-seq marks from the ENCODE 
export work=/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/data
export script=/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step2/script

bedtools intersect -a /public/home/ccc2018/data/database/Encode/Lung/download/bed/all_lung_peak.bed -b ${work}/sig05_case_AD_SC_60bp.bed -wa -wb > ${work}/sig05_case_AD_SC_60bp_EncodeLung.bed

bedtools intersect -a /public/home/ccc2018/data/database/Encode/Lung/download/bed/all_lung_peak.bed -b ${work}/sig05_case_AD_SC_vid.bed -wa -wb > ${work}/sig05_case_AD_SC_vid_EncodeLung.bed

####Annotation with eQTL in lung tissues
Rscript ${script}/eqtl_NMU.R

###Predicted TFBS disruptions
Rscript ${script}/motifbreakR.R
TFBSs disruptions were also scored using the SNP2TFBS webtool(https://epd.expasy.org/snp2tfbs)

###Annotation with LungENN
###model traning and prediction using selene
##The algorithmic framework is available at https://github.com/FunctionLab/selene/;
#The predicted results of LungENN are in the directory of ./data/LungENN/ 
Rscript ${script}/LungENN.R

###combind MPRA results and annotations
Rscript ${script}/regulatory_exp_anno.R
Rscript ${script}/regulatory_allele_anno.R
Rscript ${script}/regulatory_allele.R

###enrichment and predcited score of frVars(figure 2)
Rscript ${script}/regulatory_exp.R
Rscript ${script}/reg_enrich_exp.R

###statistic analysis of fvars and causal variants identification (data for table 1, table 2 and figure 2)
Rscript ${script}/fvars_stat.R

#####Step 3. plot
export work=/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step3/data
export script=/public/home/ccc2018/data/project/MPRA/LC/GWAS/revise/1/data_for_revision/manuscript/step3/script

####fig2a-h
Rscript ${script}/fig2.R

####fig3
###fig3a
locuszoom=/path/to/locuszoom/
locuszoom=/path/to/locuszoom/

python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/data_4q22.1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/LD_4q22.1.txt --build hg19 --chr 4 --start 89750000 --end 90100000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/bed_4q22.1.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'

python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/GWASdata_4q22.1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/LD_4q22.1.txt --build hg19 --chr 4 --start 89750000 --end 90100000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/bed_4q22.1.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'

###fig3b,3f,3g
experiments data were presented in source_data.excel

###fig3c
Rscript ${script}/motifbreakR.R

###fig3d
summary data for eqtl was presented in sig05_case_AD_SC_vid_NMU_Lung_eqtl.txt
Individual data for eqtl was presented in source_data.csv
Rscript ${script}/plot_eQTL.R

###fig3e
Rscript ${script}/plot_COLOC_edit.R

###fig3e
Rscript ${script}/fig3h.R

####fig4
locuszoom=/path/to/locuszoom/
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig4a/data_5p15.33.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig4a/LD_5p15.33.txt --build hg19 --chr 5 --start 1249000 --end 1361000 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig4a/bed_5p15.33.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'

python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/GWASdata_5p15.33.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/LD_5p15.33.txt --build hg19 --chr 5 --start 1249000 --end 1361000 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig3a/bed_5p15.33.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'

###fig.4c
Rscript ${script}/fig4c.R

###fig4d-4H
experiments data were presented in source_data.excel
summary data for eqtl was presented in sig05_case_AD_SC_vid_NMU_Lung_eqtl.txt

####fig5
locuszoom=/path/to/locuszoom/
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig5a/GWASdata_20q11.23.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig5a/LD_20q11.23.txt --build hg19 --chr 20 --start 35490000 --end 35600000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig5a/bed_20q11.23.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'
  python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig5a/data_20q11.23.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig5a/LD_20q11.23.txt --build hg19 --chr 20 --start 35490000 --end 35600000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/fig5a/bed_20q11.23.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'

###fig5b-5e
experiments data were presented in source_data.excel
summary data for eqtl was presented in sig05_case_AD_SC_vid_NMU_Lung_eqtl.txt

######fig6
Rscript ${script}/fig6.R

######figs1
Rscript $script/OligoCountsSumStat_rep6.R ./

######figs2
sh ${script}/figs2.sh
Rscript ${script}/figs2.R
Rscript ${script}/plot_TFBS_corr.R

######figs3
Rscript ${script}/plot_pred_score.R

#####figs4-12
locuszoom==/path/to/locuszoom/
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/GWASdata_3q28.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_3q28.txt --build hg19 --chr 3 --start 189210844 --end 189460844 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_3q28.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_3q28.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_3q28.txt --build hg19 --chr 3 --start 189210844 --end 189460844 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_3q28.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
##figs5
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_10q25.2.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_10q25.2.txt --build hg19 --chr 10 --start 114400000 --end 114600000 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_10q25.2.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none' &
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/GWASdata_10q25.2.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_10q25.2.txt --build hg19 --chr 10 --start 114400000 --end 114600000 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_10q25.2.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
######figs6
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/GWASdata_17q24.2.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_17q24.2.txt --build hg19 --chr 17 --start 65714382 --end 66163360 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_17q24.2.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_17q24.2.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_17q24.2.txt --build hg19 --chr 17 --start 65714382 --end 66163360 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_17q24.2.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
######figs7
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/GWASdata_14q13.1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_14q13.1.txt --build hg19 --chr 14 --start 35190000 --end 35390000 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_14q13.1.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_14q13.1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_14q13.1.txt --build hg19 --chr 14 --start 35190000 --end 35390000 --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_14q13.1.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
#####figs8-1
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/GWASdata_11q23.3_1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_11q23.3_1.txt --build hg19 --chr 11 --start 118030000 --end 118150000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_11q23.3.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_11q23.3_1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_11q23.3_1.txt --build hg19 --chr 11 --start 118030000 --end 118150000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_11q23.3.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
#####figs8-2
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/GWASdata_11q23.3_1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_11q23.3_1.txt--build hg19 --chr 11 --start 119080000 --end 119151000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_11q23.3.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_11q23.3_1.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_11q23.3_1.txt --build hg19 --chr 11 --start 119080000 --end 119151000  --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_11q23.3.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none'&
######figs9
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_8p12.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_8p12.txt --build hg19 --chr 8 --start 32205979 --end 32605979 --flank 200kb --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_8p12.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="2" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none' ymax=8&
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/GWASdata_8p12.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_8p12.txt --build hg19 --chr 8 --start 32205979 --end 32605979 --flank 200kb --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_8p12.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none' ymax=8&
#####figs10
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_15q23.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_15q23.txt --build hg19 --refsnp 15:69593622 --flank 250kb --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_15q23.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none' ymax=8&
####figs11
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_4p15.33.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_4p15.33.txt --build hg19 --refsnp 4:18022834 --flank 500kb --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_4p15.31.txt --plotonly \
--markercol SNP --pvalcol P signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none' ymax=8&
######figs12
python2.7 ${locuszoom} --metal /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/data_6p21.2.txt --ld /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/LD_6p21.2.txt --build hg19 --refsnp 6:40467616 --flank 500kb --bed-tracks /public/home/ccc2018/data/project/MPRA/LC/GWAS/check_data_20240208/data/figs4-s12/bed_6p21.2.txt --plotonly \
--markercol SNP --pvalcol P  signifLine="5" signifLineColor="darkgreen" \
--gene-table gencode --no-cleanup showRecomb=F width=12 height=10 axisTextSize=2 axisSize=2 legend='none' ymax=8&
