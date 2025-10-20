rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./19_TIDE")){
  dir.create("./19_TIDE")
}
setwd("./19_TIDE")

#Sort out the data-----------------
dat.tcga <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")
dat.tcga <- dat.tcga[,riskScore$sample]

###TIDE------------------------------------------------
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# tide_dat[which(tide_dat<0)] <- 0
tide_dat <- apply(dat.tcga,2,FPKM2TPM)
tide_dat <- log2(tide_dat + 1)
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)

####The website conducts TIDE prediction: http://tide.dfci.harvard.edu/


tide_result <- read.csv("tide_result.csv",header = T)
colnames(tide_result)
#View(tide_result)
tide_result2 <- subset(tide_result, select = c("Patient","Responder","TIDE","Dysfunction","Exclusion"))
tide_plot_dat <- data.frame(Patient=riskScore$sample,
                            riskScore=riskScore$riskScore,
                            risk_group=riskScore$risk)
tide_plot_dat <- merge(tide_plot_dat, tide_result2, by = "Patient")
tide_plot_dat$risk_group <- factor(tide_plot_dat$risk_group,
                                   levels = c("High Risk","Low Risk"))

dim(tide_plot_dat)
colnames(tide_plot_dat)
class(tide_plot_dat$riskScore)
# tide_plot_dat$riskScore <- as.numeric(tide_plot_dat$riskScore)

data2 <- tide_plot_dat[,-c(2,4)]
data2<-gather(data2,
              key = tide.type,
              value = score,
              -c("Patient",'risk_group'))
data2 <- data.frame(data2)


# control.sample <- group$sample[which(group$group=='control')]
# ## Sample grouping
# data2$Group<-ifelse(data2$sample%in%control.sample,'control','COPD')
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-data2%>%
  group_by(tide.type)%>%
  wilcox_test(score ~ risk_group)%>%
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
data2$risk_group <- factor(data2$risk_group,levels = c("High Risk","Low Risk"))
exp_plot <- ggplot(data2,aes(x = risk_group, y = score, color = risk_group)) +
  geom_violin(trim=F,color='black',aes(fill=risk_group)) + # Draw a violin picture. Use "color=" to set the color of the outline of the violin picture. (# If you don't want the outline, you can set it to white. Below, setting the background to white actually indicates that you don't want the outline.)
  #If "trim" is TRUE(the default value), the tail of the violin will be trimmed to the data range. If it is FALSE, do not trim the tail.
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9),
               show.legend = F) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,
               show.legend = F)+ #Draw the box plot. Here, width=0.1 controls the width of the box plot in the violin plot
  scale_fill_manual(values= c("#106eb4","#e3ba20"), name = "risk")+
  labs(title="TIDE", x="", y = "",size=20) +
  scale_colour_manual(values = c("#106eb4","#e3ba20"))+
  # scale_color_brewer(palette = 'Set2')+
  # stat_compare_means(aes(group = risk,label = "p.format"),
  #                    method = "wilcox.test",label = 'p.signif',label.x = 1.4)+
  stat_pvalue_manual(stat.test,
                     y.position = 3.5,
                     size = 3.2,
                     family = "Times",
                     label = "p1",  #Select the results from the difference analysis or the W-test
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


#Statistical bar chart --------

# library
data3 <- tide_plot_dat[,c(1,3,4)]


data3 <- data.frame(sample=data3$Patient, risk=data3$risk_group, Responder=data3$Responder)
rep1 <- xtabs(~risk+Responder,data=data3)
p <- chisq.test(rep1)[[3]]

###draw
rep2 <- table(data3)%>%as.data.frame()
## Stacked graph
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
  labs(x="",y="Percentage",fill="Responder",title = "Responder*")
p1
ggsave(filename = '02.TIDE_barchart.pdf',p1,w=6,h=5)
ggsave(filename = '02.TIDE_barchart.png',p1,w=6,h=5,dpi = 600)



###IPS--------------------------------------
ips<-read_delim('TCIA-ClinicalData.tsv')%>%data.frame(.)
ips<-ips[,c('barcode', 'ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos', 'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')]
rownames(ips)<-ips$barcode

risk <- riskScore
risk$sample<-substr(risk$sample,1,12)

risk<-risk[!duplicated(risk$sample),] #Keep the non-repetitive ones because they are all tumor samples and generally there won't be any duplicates
rownames(risk)<-risk$sample
data<-merge(risk,ips,by = 'row.names')
rownames(data)<-data$Row.names
# write.csv(data,'02.ips_data.csv',quote=F) 
#Difference analysis------------
library(ggplot2)
library(ggpubr)
#Draw all four fractions
plot.ips <- c('ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos', 'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')

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
  
  png(paste0('03.',plot.ips[i],'.png'),width=7,height=6,family='Times',units='in',res=600)
  print(p)
  dev.off()
  pdf(paste0('03.',plot.ips[i],'.pdf'),width=7,height=6,family='Times')
  print(p)
  dev.off()
}

