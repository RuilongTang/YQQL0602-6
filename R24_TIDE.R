rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./19_TIDE")){
  dir.create("./19_TIDE")
}
setwd("./19_TIDE/new_TIDE/")

#整理数据-----------------
# dat.tcga <- read.csv("../../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
# riskScore <- readRDS("../../08_LASSO/riskScore.rds")
# riskScore$sample <- rownames(riskScore)
# riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")
# dat.tcga <- dat.tcga[,riskScore$sample]
# 
# ###TIDE------------------------------------------------
# FPKM2TPM <- function(fpkm){
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# # tide_dat[which(tide_dat<0)] <- 0
# tide_dat <- apply(dat.tcga,2,FPKM2TPM)
# tide_dat <- log2(tide_dat + 1)
# rownmean <- apply(tide_dat,1,mean)
# tide_dat2 <- sweep(tide_dat, 1, rownmean)
# dim(tide_dat2)
# write.table(tide_dat2,
#             file ="tide_dat.txt",
#             sep = "\t",
#             quote = F,
#             row.names = T)
# 
# ####网站TIDE预测http://tide.dfci.harvard.edu/
# 


msi <- read.csv("../../18_MSI/MSI.csv",header = T)
tide_result <- read.csv("../tide_result.csv",header = T)
colnames(tide_result)
#View(tide_result)
tide_result2 <- subset(tide_result, select = c("Patient","Responder","TIDE","Dysfunction","Exclusion"))
rownames(tide_result2) <- tide_result2$Patient
tide_result2 <- tide_result2[msi$sample,]

tide_plot_dat <- data.frame(Patient=tide_result2$Patient,
                            riskScore=msi$riskScore,
                            risk_group=msi$risk,
                            group = msi$group)
tide_plot_dat <- merge(tide_plot_dat, tide_result2, by = "Patient")
tide_plot_dat$group <- factor(tide_plot_dat$group,
                                   levels = c("MSI-H","MSI-L/MSS"))

dim(tide_plot_dat)
colnames(tide_plot_dat)
class(tide_plot_dat$riskScore)
# tide_plot_dat$riskScore <- as.numeric(tide_plot_dat$riskScore)

data2 <- tide_plot_dat[,-c(2,3,5)]
data2<-gather(data2,
              key = tide.type,
              value = score,
              -c("Patient",'group'))
data2 <- data.frame(data2)


# control.sample <- group$sample[which(group$group=='control')]
# ## 样本分组
# data2$Group<-ifelse(data2$sample%in%control.sample,'control','COPD')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-data2%>%
  group_by(tide.type)%>%
  wilcox_test(score ~ group)%>%
  adjust_pvalue(method = 'fdr')

# diff <- read.csv('../01_DEGs/01.DEG_sig.csv',row.names = 1)
# df<-diff[stat.test$Symbol,]
# stat.test$change<-df$change#%>%as.numeric()%>%round(digits = 3)
stat.test$p1<-ifelse(stat.test$p<0.0001,"****",
                     ifelse(stat.test$p<0.001,"***",
                            ifelse(stat.test$p<0.01,"**",
                                   ifelse(stat.test$p<0.05,"*",'ns'))))
library(ggplot2)
library(ggpubr)
data2$group <- factor(data2$group,levels = c("MSI-H","MSI-L/MSS"))
exp_plot <- ggplot(data2,aes(x = group, y = score, color = group)) +
  geom_violin(trim=F,color='black',aes(fill=group)) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9),
               show.legend = F) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,
               show.legend = F)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#106eb4","#e3ba20"), name = "Group")+
  labs(title="TIDE", x="", y = "",size=20) +
  scale_colour_manual(values = c("#106eb4","#e3ba20"))+
  # scale_color_brewer(palette = 'Set2')+
  # stat_compare_means(aes(group = risk,label = "p.format"),
  #                    method = "wilcox.test",label = 'p.signif',label.x = 1.4)+
  stat_pvalue_manual(stat.test,
                     y.position = 3.5,
                     size = 3.2,
                     family = "Times",
                     label = "p1",  #选择差异分析中的结果或者w检验
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family='Times'))+
  facet_wrap(~tide.type,scales = "free",nrow = 1) 
exp_plot
ggsave(filename = '01.TIDE.pdf',exp_plot,w=8,h=4)
ggsave(filename = '01.TIDE.png',exp_plot,w=8,h=4,dpi = 600)


#统计条形图--------

# library
data3 <- tide_plot_dat[,c(1,4,5)]


data3 <- data.frame(sample=data3$Patient, risk=data3$group, Responder=data3$Responder)
rep1 <- xtabs(~risk+Responder,data=data3)
p <- chisq.test(rep1)[[3]]

###绘图
rep2 <- table(data3)%>%as.data.frame()
## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)

p1 <- ggplot(rep2,aes(x=risk,y=Freq,fill=Responder))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(3, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=30, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="Responder",title = "Responder")
p1
ggsave(filename = '02.TIDE_barchart.pdf',p1,w=6,h=5)
ggsave(filename = '02.TIDE_barchart.png',p1,w=6,h=5,dpi = 600)



###IPS MSS--------------------------------------
ips<-read_delim('../TCIA-ClinicalData.tsv')%>%data.frame(.)
ips<-ips[,c('barcode', 'ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos', 'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')]
rownames(ips)<-ips$barcode

risk <- msi
risk <- risk[risk$group == "MSI-L/MSS",]
risk$sample<-substr(risk$sample,1,12)

risk<-risk[!duplicated(risk$sample),] #保留不重复的，因为都是肿瘤样本，一般不会有重复的
rownames(risk)<-risk$sample
data<-merge(risk,ips,by = 'row.names')
rownames(data)<-data$Row.names
# write.csv(data,'02.ips_data.csv',quote=F) 
#差异------------
library(ggplot2)
library(ggpubr)
#4个分数都画
plot.ips <- c('ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos', 'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')


setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/19_TIDE/new_TIDE/MSS-risk/")
for (i in c(1:4)) {
  # i <- 1
  my_comparisons = list( c("High Risk","Low Risk"))
  p <- ggplot(data, aes(risk, get(plot.ips[i]))) + 
    ggdist::stat_halfeye(aes(color=risk,fill=risk),adjust = .5, width = .7, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(aes(color=risk),width = .1, outlier.shape = NA) +
    # gghalves::geom_half_point(aes(color=Species),side = "l", range_scale = .4, alpha = .5) +
    #ggsci::scale_color_nejm()+
    #ggsci::scale_fill_nejm() +
    geom_jitter(aes(color=risk),width = .05, alpha = .3) +
    scale_fill_manual(values=c("#e57667","#23b3bb")) +
    scale_color_manual(values=c("#e57667","#23b3bb")) + 
    coord_flip()+
    theme_bw()+
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "wilcox.test",label.y = c(10.5),
                       color='black',family="Times",face = "bold")+
    labs(x = "risk", y = plot.ips[i], title = "") + 
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times",color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
          axis.text.y = element_text(size = 0, face = "bold", family = "Times",color='black'),
          legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
          legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
          legend.position = "bottom"
    )+
    theme(panel.grid =element_blank()) +   
    theme(axis.text = element_blank()) +  
    theme(axis.ticks = element_blank())+
    theme(panel.background = element_blank())
  
  png(paste0(plot.ips[i],'.png'),width=7,height=6,family='Times',units='in',res=600)
  print(p)
  dev.off()
  pdf(paste0(plot.ips[i],'.pdf'),width=7,height=6,family='Times')
  print(p)
  dev.off()
}



###IPS MSS-MSI--------------------------------------
risk <- msi
risk$sample<-substr(risk$sample,1,12)

risk<-risk[!duplicated(risk$sample),] #保留不重复的，因为都是肿瘤样本，一般不会有重复的
rownames(risk)<-risk$sample
data<-merge(risk,ips,by = 'row.names')
rownames(data)<-data$Row.names
# write.csv(data,'02.ips_data.csv',quote=F) 
#差异------------
library(ggplot2)
library(ggpubr)
#4个分数都画
plot.ips <- c('ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos', 'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/19_TIDE/new_TIDE/MSI/")
for (i in c(1:4)) {
  # i <- 1
  my_comparisons = list( c("MSI-H","MSI-L/MSS"))
  p <- ggplot(data, aes(group, get(plot.ips[i]))) + 
    ggdist::stat_halfeye(aes(color=group,fill=group),adjust = .5, width = .7, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(aes(color=group),width = .1, outlier.shape = NA) +
    # gghalves::geom_half_point(aes(color=Species),side = "l", range_scale = .4, alpha = .5) +
    #ggsci::scale_color_nejm()+
    #ggsci::scale_fill_nejm() +
    geom_jitter(aes(color=group),width = .05, alpha = .3) +
    scale_fill_manual(values=c("#e57667","#23b3bb")) +
    scale_color_manual(values=c("#e57667","#23b3bb")) + 
    coord_flip()+
    theme_bw()+
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "wilcox.test",label.y = c(10.5),
                       color='black',family="Times",face = "bold")+
    labs(x = "risk", y = plot.ips[i], title = "") + 
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times",color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
          axis.text.y = element_text(size = 0, face = "bold", family = "Times",color='black'),
          legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
          legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
          legend.position = "bottom"
    )+
    theme(panel.grid =element_blank()) +   
    theme(axis.text = element_blank()) +  
    theme(axis.ticks = element_blank())+
    theme(panel.background = element_blank())
  
  png(paste0(plot.ips[i],'.png'),width=7,height=6,family='Times',units='in',res=600)
  print(p)
  dev.off()
  pdf(paste0(plot.ips[i],'.pdf'),width=7,height=6,family='Times')
  print(p)
  dev.off()
}

