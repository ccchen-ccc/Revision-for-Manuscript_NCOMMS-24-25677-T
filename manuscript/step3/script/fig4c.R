####fig4c
library(forestplot)
labeltext <-  cbind(c("statu",'LUAD','LUSC','no-smoker','smoker'),
                    c("Odds Ratio (95% CI)",'1.39(1.34-1.45)','1.06(1.01-1.13)','1.39(1.32-1.46)','1.21(1,15-1.27)'))
mean=c(1.39,1.06,1.39,1.21)
lower=c(1.34,1.01,1.32,1.15)
upper=c(1.45,1.13,1.46,1.27)
pdf(file = './fig4c',width = 10,height = 8)
forestplot(labeltext=labeltext,
           graph.pos=2, #为Pvalue箱线图所在的位置
           mean=c(NA,mean),
           lower=c(NA,lower), upper=c(NA,upper),
           #定义标题
           title="Odds Ratio Plot",
           ##定义x轴
           #xlab="",
           ##根据亚组的位置，设置线型，宽度造成“区块感”
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           ##fpColors函数设置颜色
           col=fpColors(box=c("#ff7d00"), lines="#ff7d00", zero = "gray50"),
           #箱线图中基准线的位置
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱子大小，线的宽度
           lwd.ci=3, boxsize=0.15
           #箱线图两端添加小竖线，高度
           #ci.vertices=TRUE, ci.vertices.height = 0.4
)
dev.off()#######输出后ai美化