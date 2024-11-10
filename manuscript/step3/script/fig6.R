####fig6
library(forestplot)
library(tidyverse)
library(readxl)
forest_plot1 <- readxl::read_xlsx('PRS_1216.xlsx',sheet = 1)
forest_plot1 <- forest_plot1 %>% 
  separate(PRS_inter_OR_CI_cox,into = c('OR','CI'),sep = ' ',remove = F) %>% 
  mutate(CI=str_replace(CI,pattern = '\\(',replacement = '')) %>% 
  mutate(CI=str_replace(CI,pattern = '\\)',replacement = '')) %>% 
  separate(CI,into = c('low','up'),'-')
labeltext <-  cbind(c("type",forest_plot1$type),
                    c("Odds Ratio (95% CI)",forest_plot1$PRS_inter_OR_CI_cox),
                    c("p",signif(forest_plot1$PRS_inter_P_cox,digits = 2)))
pdf('./fig6.pdf',width = 10,height = 8)
forestplot(labeltext=labeltext,
           graph.pos=2, 
           mean=c(NA,as.numeric(forest_plot1$OR)),
           lower=c(NA,as.numeric(forest_plot1$low)), upper=c(NA,as.numeric(forest_plot1$up)),
           title="Odds Ratio Plot",
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           col=fpColors(box=c("#ff7d00"), lines="#ff7d00", zero = "gray50"),
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           lwd.ci=3, boxsize=0.15)
dev.off()
