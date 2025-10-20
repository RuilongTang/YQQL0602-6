rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./22_Expression")){
  dir.create("./22_Expression")
}
setwd("./22_Expression")


## 01 Read data ----------
data <- read.csv("../00_Rowdata/fpkm.tcga.csv",header=TRUE, row.names=1, check.names=FALSE)
data <- log2(data+1)
group <- read.csv("../00_Rowdata/tcga.group.csv",header = T,row.names = 1)
control <- group$id[group$group == "Control"]
gene <- readRDS("../08_LASSO/cox_gene.rds")
gene_fpkm <- data[gene,]

gene_fpkm$gene<-rownames(gene_fpkm)
violin_dat <- gather(gene_fpkm, key=sample, value='expr', -c("gene"))
violin_dat$group <- ifelse(violin_dat$sample %in% control,"Control","CRC")

violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=expr,
                                      fill=group)) +
  #geom_violin(trim=F,color="black") + # Draw a violin picture. Use "color=" to set the color of the outline of the violin picture. (# If you don't want the outline, you can set it to white. Below, setting the background to white actually indicates that you don't want the outline.)
  #If "trim" is TRUE(the default value), the tail of the violin will be trimmed to the data range. If it is FALSE, do not trim the tail.
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               aes(fill=group),
               position=position_dodge(0.9),
               outlier.shape = NA)+ #Draw the box plot. Here, width=0.1 controls the width of the box plot in the violin plot
  scale_fill_manual(values= c("#658fcb", "#f38687"), name = "Group")+
  labs(title="Expression level", x="", y = "Log2(Exp+1)",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  # stat_pvalue_manual(DE.checkpoint,
  #                    x ="gene",
  #                    y.position = 7,
  #                    size = 3.2,
  #                    color = "black",
  #                    family = "Times",
  #                    label = "p.adj")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 14),
        text = element_text(family = 'Times'),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#+
#facet_wrap(~gene,scales = "free",nrow = 2) #Separate or combine them together, two lines or one line
violin_plot

ggsave('01.expression.pdf',violin_plot,w=10,h=6)
ggsave('01.expression.png',violin_plot,w=10,h=6)






