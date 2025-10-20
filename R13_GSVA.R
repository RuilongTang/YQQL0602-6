rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./13_GSVA")){
  dir.create("./13_GSVA")
}
setwd("./13_GSVA")



## 01 Read data----------
data <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")
data <- log2(data+1)
dat <- data[,riskScore$sample]

#Background gene set
library(GSVA)
library(GSEABase)
library(limma)
KEGG_ref <- getGmt("../../../pipeline/GSVA/h.all.v2023.1.Hs.symbols.gmt")

#GSVA analysis
es_KEGG <- gsva(as.matrix(dat), KEGG_ref,
                min.sz=2, max.sz=1000, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)

saveRDS(es_KEGG,"es_KEGG.rds")
es_KEGG <- readRDS("es_KEGG.rds")




#Expression and grouping files----------------
group <- riskScore[,c('sample','risk')]%>%data.frame()
colnames(group) <- c('id','group')
group <- group[order(group$group,decreasing = TRUE),]%>%data.frame()
design <- model.matrix(~0 + group$group)
rownames(design) <- group$id
colnames(design) <- c("High","Low")
compare<- makeContrasts("High-Low", levels = design)


#Difference analysis--------------------------
fit <- lmFit(es_KEGG, design)
fit2 <- contrasts.fit(fit ,compare)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

#Differential results
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$t) > logFCcutoff,
         ifelse(allGeneSets$t > logFCcutoff,'UP','DOWN'),'NoSignifi')
)
write.csv(allGeneSets, file ="01.GSVA_Hallmark_all.csv",row.names = T)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$t) > logFCcutoff )
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# 28  7
write.csv(DEGeneSets,file = "01.GSVA_Hallmark_sig.CSV",row.names = T)

##Heat map------------------
library(pheatmap)

group.gsva<-group
DEGeneSets <- DEGeneSets[order(DEGeneSets$t,decreasing = T),]
head_5 <- head(DEGeneSets,5) %>% rownames(.)
tail_5 <- tail(DEGeneSets,5) %>% rownames(.)
kegg_path <- c(head_5,tail_5)
heat<-es_KEGG[kegg_path,] ## Score according to the groups


annotation_col <- as.data.frame(group.gsva$group)
colnames(annotation_col) = 'Group'
rownames(annotation_col) = group.gsva$id
heat <- heat[,rownames(annotation_col)]

color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(Group = c('High Risk'="#FFAEB9",'Low Risk'="#00FFFF"))

pdf(file = "02.heatmap.pdf",w=11,h=7,family = 'Times')
pheatmap::pheatmap(heat,
                   color = colorRampPalette(color.key)(50),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellheight=40,
                   labels_row = NULL,
                   clustering_method = 'ward.D2',
                   show_rownames = T,
                   show_colnames = F,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = T,
                   width = 10,
                   main = "GSVA Hallmark")
dev.off()

png(file = "02.heatmap.png",w=11,h=7,units = 'in',res = 300, family = 'Times')
pheatmap::pheatmap(heat,
                   color = colorRampPalette(color.key)(50),
                   border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellheight=40,
                   labels_row = NULL,
                   clustering_method = 'ward.D2',
                   show_rownames = T,
                   show_colnames = F,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = T,
                   width = 10,
                   main = "GSVA Hallmark")
dev.off()

### Divergent bar chart drawing (Differences)-------------

library(ggprism)
## barplot
dat_plot<-data.frame(id=rownames(DEGeneSets),
                     t=DEGeneSets$t)
DEGeneSets <- DEGeneSets[dat_plot$id,]
dat_plot$threshold = factor(DEGeneSets$change)
# dat_plot$threshold = factor(ifelse(dat_plot$p  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
table(dat_plot$threshold)
dat_plot$threshold <- ifelse(dat_plot$threshold=='UP','Up',
                             ifelse(dat_plot$threshold=='DOWN','Down','NoSignifi'))
dat_plot<-dat_plot%>%arrange(t)

dat_plot$id<-factor(dat_plot$id,levels=dat_plot$id)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#a64545','NoSignifi'='#cccccc','Down'='#779fc2')) +
  #  geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score') + #Note that the coordinate axes have rotated
  ##guides(fill="none")+ # No legend is displayed
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
low1 <- dat_plot %>% filter(t < -0) %>% nrow()
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 0) %>% nrow()
high1 <- nrow(dat_plot)

# Add labels from bottom to top in sequence
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black',family='Times')  + # Labels less than -1 are black
  # geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
  #            hjust = 0,color = 'grey') + # Gray label
  # geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
  #            hjust = 1,color = 'grey') + # Gray label
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black',family='Times') # Labels more than -1 are black
p

png(file = "03.Bar.png",w=10,h=6,units = 'in',res = 300, family = 'Times')
p
dev.off()
pdf(file = "03.Bar.pdf",w=10,h=6, family = 'Times')
p
dev.off()

