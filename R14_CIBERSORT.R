rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./14_CIBERSORT")){
  dir.create("./14_CIBERSORT")
}
setwd("./14_CIBERSORT")

## 01 Read data----------
data <- read.csv("../00_Rowdata/tumor.tcga.csv",header=TRUE, row.names=1, check.names=FALSE)
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")
data <- log2(data+1)
dat <- data[,riskScore$sample]


####Immune-infiltrating cells --------------------------------------------
library(immunedeconv)
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")
res.cibersort<-deconvolute(as.matrix(dat),method = 'cibersort')
write.csv(res.cibersort,'cibersort.csv')
saveRDS(res.cibersort,"cibersort.rds")
res.cibersort <- readRDS("cibersort.rds")


####Draw the petal diagram------------------------------------------
library(ggplot2)
group <- riskScore[,c(2:3)]
high <- group$sample[which(group$risk =="High Risk")]
low <- group$sample[which(group$risk =="Low Risk")]
dat_ciber <- res.cibersort
dat_ciber$high <- rowSums(res.cibersort[high])
dat_ciber$low <- rowSums(res.cibersort[low])
result <- dat_ciber[,c(1,608,609)]  ###Extract the columns "high" and "low"
###Delete rows that are all zeros
result <- result[which(rowSums(result==0)==0),]


# Get the name of each lable and its position on the Y-axis
label_data <- result
label_data$id <- seq(1, nrow(label_data))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 0, 0)
# label_data$angle <- ifelse(angle < -90, angle+180, angle)
label_data$angle <- ifelse(angle < -90, angle-70, angle-70)

###petal diagram---------------
df <- result
# df <- df[order(-df$high),]
x <- 1:(21*18)  ###Sample size*18
y<-sin(21*x*pi/(21*18))
df1<-data.frame(x1=x,
                y1=abs(y),
                var=gl(21,18,labels = df$cell_type))

merge(df1,df,by.x = 'var',by.y = 'cell_type') %>% 
  mutate(new_high=y1*high,new_low=y1*low) -> df2

df2$cell_type <- factor(df2$var,levels = df$cell_type)

###high
p1 <- ggplot(data=df2,aes(x=x1,y=new_high))+
  geom_area(aes(fill=cell_type),
            alpha=0.8,
            color="black",
            show.legend = T)+
  coord_polar()+
  geom_text(data=label_data, aes(x=id*18-18, y=70, label=cell_type, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=4, 
            angle= label_data$angle, inherit.aes = FALSE )+
  theme_bw()+
  ylim(0,100)+
  theme(axis.text.x = element_text(size = 0,face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=10),
        legend.title = element_text(face = "bold", hjust = 0.5,colour="black", size=15),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = c("white","grey","white","grey","white","grey")))
p1
ggsave(filename = '01.high.pdf',p1,w=14,h=11)
ggsave(filename = '01.high.png',p1,w=14,h=11,bg = "white") 

###low
p1 <- ggplot(data=df2,aes(x=x1,y=new_low))+
  geom_area(aes(fill=cell_type),
            alpha=0.8,
            color="black",
            show.legend = T)+
  coord_polar()+
  geom_text(data=label_data, aes(x=id*18-18, y=70, label=cell_type, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=4, 
            angle= label_data$angle, inherit.aes = FALSE )+
  theme_bw()+
  ylim(0,100)+
  theme(axis.text.x = element_text(size = 0,face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=10),
        legend.title = element_text(face = "bold", hjust = 0.5,colour="black", size=15),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = c("white","grey","white","grey","white","grey")))
p1
ggsave(filename = '02.low.pdf',p1,w=14,h=11)
ggsave(filename = '02.low.png',p1,w=14,h=11,bg = "white") 



## Difference results ---------------
dat.cibersort <- res.cibersort %>% 
  tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
group <- riskScore[,c(2,3)]

dat.cibersort <- merge(group, dat.cibersort, by = "sample")
###Delete all columns whose values are all 0 purrr::discard(df, ~all(.x == 0)) purrr::keep(df, ~any(.x != 0))
dat.cibersort <- dat.cibersort[,apply(dat.cibersort,2,function(dat.cibersort) !all(dat.cibersort==0))]
dat.cibersort2 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "risk"))
library(rstatix)
library(ggplot2)
library(ggpubr)
colnames(dat.cibersort2)
stat_cibersort <- dat.cibersort2 %>% 
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ risk) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.csv(stat_cibersort,file = 'stat.cibersort.csv')
DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.csv(DE.cibersort,file = 'DE.cibersort.csv')
colnames(dat.cibersort2)

# violin.cibersort<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]
violin.cibersort<-dat.cibersort2

violin_plot <- ggplot(violin.cibersort, aes(x=ImmuneCell, y=Score,fill=risk)) +
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  # geom_violin(aes(fill=group),trim=FALSE,color="black") + #To draw a violin picture, use "color=" to set the color of the outline of the violin picture. (# No outline can be set to white. Below, setting the background to white actually indicates no outline.)
  #If "trim" is TRUE(the default value), the tail of the violin will be trimmed to the data range. If it is FALSE, do not trim the tail.
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #Draw the box plot. Here, width=0.1 controls the width of the box plot in the violin plot
  # geom_point(aes(color = group), position=position_dodge(0.9), size = 1) +
  scale_fill_manual(values= c("#DC143C","#1E90FF"), name = "Group")+
  labs(title="Immune Cell", x="", y = "",size=20) +
  stat_compare_means(data = violin.cibersort,
                     mapping = aes(group = risk),
                     label ="p.signif",
                     method = 'wilcox.test',
                     label.y = 0.7,
                     size = 6,
                     paired = F) +
  theme_bw()+
  labs(y = "Immune score")+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot
ggsave(filename = '04.violin.plot.pdf',violin_plot,w=12,h=8)
ggsave(filename = '04.violin.plot.png',violin_plot,w=14,h=8)
dev.off()



###The correlation between risk score and immune cells --------------------------------------------
library(psych)


tiic <- res.cibersort%>%column_to_rownames(var = 'cell_type')
tiic <- tiic[DE.cibersort$ImmuneCell,]%>%t%>%as.data.frame()

x <-as.numeric(riskScore[,1])
y <-tiic
d <- corr.test(y,x,use="complete",method = 'spearman')
r <- data.frame(d$r)
p <- data.frame(d$p)
correlation<-data.frame(rownames(p),r$d.r,p$d.p)
colnames(correlation) <- c("cell","Correlation","p.value")
correlation <- na.omit(correlation)
correlation<-correlation[correlation$p.value<0.05,]
correlation$sig[correlation$p.value>0.05] <- "ns"   
correlation$sig[correlation$p.value<0.05&correlation$p.value>0.01] <- "*"   
correlation$sig[correlation$p.value<0.01&correlation$p.value>0.001] <- "**"  
correlation$sig[correlation$p.value<0.001&correlation$p.value>0.0001] <- "***"   
correlation$sig[correlation$p.value<0.0001] <- "****"
correlation$'correlation'<-abs(correlation$Correlation)
library("viridis")
#trace(ggdotchart,edit=T)
p0<-ggdotchart(correlation, x = "cell", y = "Correlation",
               dot.size ='correlation',
               color ='p.value',
               sorting = "descending",
               add = "segments",                             # Add a stick
               rotate = TRUE,
               ggtheme = theme_pubr(),                        # Change the topic
               ylab="Correlation Coefficient (r)",
               xlab='',
               # dot.size = 6 ,
               title="RiskScore"
               
)+scale_colour_gradient( high = "#4682B4",low = "#66CDAA")

p10<-p0+theme(legend.position = "right",
              panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
  theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
        axis.text.x =element_text(size=12,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=12,family = "Times", face = "bold"),
        plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
        legend.text = element_text(size = 16, family = "Times"),
        legend.title = element_text(size = 18, family = "Times",face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p10
ggsave("05.Riskscore.pdf",p10,w=6,h=4)
ggsave("05.Riskscore.png",p10,w=6,h=4)


