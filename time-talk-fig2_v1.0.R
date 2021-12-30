####Fig1
####20211003
####count Ligand receptor gene ratio

####----0.1 load package--------
library(EnrichedHeatmap)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(clusterProfiler)
library(future.apply)
library(circlize)
library(ggthemes)
library(ggsci)
library(ggsignif)
library(preprocessCore)
library(factoextra)
library(Seurat)
source(file = "code/myUtils.R")


####----0.2 define color------
rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))) 
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
col.spectral <- colorRampPalette(brewer.pal(11,'Spectral')[-6])
test.color.3 <- colorRampPalette(c("#f86e11","#e9970a","#71a701","#62b474","#71c3ac","#9fc4ca"))
rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))
skyblueYelow <- colorRampPalette(c("#6db0d2","black","#f1e63d"))
skybluered <- colorRampPalette(c("skyblue","black","orange"))
solarExtra <- colorRampPalette(c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
hic.red <- colorRampPalette(c("white", "red"))
hic.orred <- colorRampPalette(colors = c(brewer.pal(9,"Reds")))
blues <- colorRampPalette(colors = brewer.pal(9,"Blues"))
ylord <- colorRampPalette(colors = brewer.pal(9,"YlOrRd"))

redwhgreen <- colorRampPalette(rev(c("#8a1719", "white", "#146533")))
mycolor.bar(redwhgreen(100),min = -1)

scales::show_col(blues(10))
###hic.pca.red <- colorRampPalette(c("#1a469c","black","#e7141a"))
hic.pca.red <- colorRampPalette(c("blue","gray1","red"))
hic.pca.redwhite <- colorRampPalette(c("#1d1856","navyblue","white","red4","#861617"))
hic.pca.orange <- colorRampPalette(c("#2f2583","black","#f9b232"))
hic.pca.skyblue <- colorRampPalette(c("skyblue","black","orange"))

#mycolor.bar(hic.pca.red(100),min = -1,max = 1)
divergentcolor <- function (n) {
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", 
                  "#984EA3", "#F29403", "#F781BF", "#BC9DCC", 
                  "#A65628", "#54B0E4", "#222F75", "#1B9E77", 
                  "#B2DF8A", "#E3BE00", "#FB9A99", "#E7298A", 
                  "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", 
                  "#B3DE69", "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

navy <- colorRampPalette(colors = c("#3e56a6", "#fffbfb", "#ee252a"))
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
coldwarm <- colorRampPalette(colors = mypalette)

bkcolor <- c(colorRampPalette(c(brewer.pal(9,"Blues" )[4:9],"#1a1919"))(50),
             colorRampPalette(c("#1a1919",rev(brewer.pal( 9,"YlOrBr" ))[1:6]))(50))

bkcolor <- colorRampPalette(colors = bkcolor)
mycolor.bar(bkcolor(10),min = -1)

coffeshop <- colorRampPalette(colors =c("#dd9b53","#d25b12","seagreen3","#c2a38a","#b8a896","grey"))
scales::show_col(colours = coffeshop(8))
scales::show_col(divergentcolor(15))

shaman.score_pal_0_col = 'lightgrey'
shaman.score_pal_pos_col = c('lightgrey', 'pink', 'red', 'orange', 'yellow')
shaman.score_pal_neg_col = c('lightgrey', 'lightblue', 'blue',"#a63ee3",'#651790')

shaman.score_pal_pos_breaks = c(1,35,50,80, 99)
shaman.score_pal_neg_breaks = c(1,35,50,80, 99)

.shaman_palette_breaks = function(n = 100, colors = shaman., breaks)
{
  colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
  for(i in 2:(length(colors)) ){
    colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
  }
  colspec = c( colspec,
               colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)])
  )
  return(colspec)
}


shaman_score_pal <- function(n=100) {
  col.neg <- rev(.shaman_palette_breaks(100,shaman.score_pal_neg_col,shaman.score_pal_neg_breaks ))
  col.pos <- .shaman_palette_breaks(100 , shaman.score_pal_pos_col,shaman.score_pal_pos_breaks )
  col.scores = colorRampPalette(colors = c(col.neg, 'lightgrey', col.pos))(n)
  
  return(col.scores)
}


mycolor.bar(my.colors = shaman_score_pal(n = 1000),min = -100,max = 100,vertical = T)
# tmp.color <- colorRampPalette(colors = c(rev(shaman.score_pal_neg_col),shaman.score_pal_0_col,shaman.score_pal_pos_col))
# mycolor.bar(my.colors = tmp.color(100),min = -100)

# scales::show_col(colours = shaman.score_pal_pos_col)
# scales::show_col(colours = shaman.score_pal_neg_col)

####-------1. QC data-----------
####------2.QC----------
####-------2.1 load data-----------

####-----2.1.1 load deng et.al data and normalized, get mean expression-----
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)
Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","fibroblast")
### test if normalized
# head(beforeEPI.meta)
# boxplot(beforeEPI.exp[,110:100])
beforeEPI.exp.mean <- beforeEPI.exp %>%
  rownames_to_column("Gene") %>%
  gather(key = "cell_id",value = "gene.exp",-Gene) %>%
  left_join(beforeEPI.meta,"cell_id") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp)) %>%
  ungroup() %>%
  mutate(Stage=factor(Stage,levels = Stage.level)) %>%
  spread(key="Stage",value = "gene.exp.mean") %>%
  column_to_rownames("Gene")

### quantile normalize
tmp.names.1 <- colnames(beforeEPI.exp.mean)
tmp.mat.1 <- preprocessCore::normalize.quantiles(as.matrix(beforeEPI.exp.mean))
colnames(tmp.mat.1) <- tmp.names.1

boxplot(tmp.mat.1)





####------2.1.2 load xiewei data--------
##### RNA-seq xiewei data
xiewei.data.exp <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myNormalize(myRemoveNA(xiewei.data.exp))


####-------2.2 scatter plot--------

####-------2.2.1 prepare for plot-----

beforeEPI.exp.mean <- readRDS(file = "res/R/beforeEPI_exp_mean_20210114.rds")
gene_symbols <- intersect(rownames(xiewei.data.exp),rownames(beforeEPI.exp.mean))

####rename for plot
colnames(beforeEPI.exp.mean)[1] <- "MII_oocyte"
colnames(xiewei.data.exp)[4] <- "late2cell"
colnames(xiewei.data.exp)[3] <- "early2cell"

print(colnames(xiewei.data.exp),quote=T)



####-------revised version-----------
Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")
data.plot.1 <- beforeEPI.exp.mean[gene_symbols,]
data.plot.2 <- xiewei.data.exp[gene_symbols,]


plot.list <- list()
for(ii in 1:length(Stage)){
  for(jj in 1:length(Stage)){
    df.plot <- data.frame(x=data.plot.1[,Stage[jj]],
                          y=data.plot.2[,Stage[ii]],
                          stringsAsFactors = F)
    res.cor <- cor(df.plot$x,df.plot$y)
    res.p <- cor.test(df.plot$x,df.plot$y)
    p <- ggplot(df.plot,aes(x,y))+
      geom_point(alpha=0.3)+
      xlab(paste0("pseudo-bulk ",Stage[jj]))+
      ylab(paste0("bulk ",Stage[ii]))+
      ggtitle(paste0("R=",round(res.cor,3)))+
      scale_x_continuous(limits = c(0,25))+
      scale_y_continuous(limits = c(0,25))+
      coord_fixed()+
      theme_cowplot()+
      theme(plot.title = element_text(hjust = 0.5))
    plot.list <- c(plot.list,list(p))
  }
}

####Figure S2
pp  <- plot_grid(plotlist = plot.list,nrow = length(Stage))
pp
myggsave(pp,prefix = "res/fig/figS2_QC_beforeEPI_cor_revised",
         suffix = ".jpg",
         width = 12,height = 12,dpi=350)


#####-----pheatmap-----

Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")
data.plot.1 <- beforeEPI.exp.mean[gene_symbols,]
data.plot.2 <- xiewei.data.exp[gene_symbols,]


res.mat <- matrix(0,nrow = length(Stage),ncol = length(Stage))
for(ii in 1:length(Stage)){
  for(jj in 1:length(Stage)){
    df.plot <- data.frame(x=data.plot.1[,Stage[jj]],
                          y=data.plot.2[,Stage[ii]],
                          stringsAsFactors = F)
    res.cor <- cor(df.plot$x,df.plot$y)
    res.p <- cor.test(df.plot$x,df.plot$y)
    res.mat[ii,jj] <- res.cor
  }
}
rownames(res.mat) <- paste0("bulk ",Stage)
colnames(res.mat) <- paste0("pseudo-bulk ",Stage)

tmp.stage  <- c("MII oocyte","zygote","early 2-cell","late 2-cell","4-cell","8-cell")

jpeg(filename = myFileName(prefix = "res/fig/Fig1_QC_pheatmap",suffix = ".jpg"),
     width = 8,height = 8,res = 350,units = "in")
ComplexHeatmap::pheatmap(res.mat,fontsize = 20,
         show_rownames = T,
         display_numbers = T,
         number_color = "white",
         fontsize_number = 24,name = "PCC",
         labels_row = paste0("bulk ",tmp.stage),
         labels_col = paste0("pseudo-bulk ",tmp.stage),
         show_colnames = T,
         heatmap_legend_param = list(title_gp = gpar(fontsize = 20,fontface = "bold"),
                                     labels_gp=gpar(fontsize = 20)),
         angle_col = "315",
         color = rdbu(100))
dev.off()

####------2.3 pca---------

####------2.3.1 replicate fig1A deng.et.al-----------

median_center <- function(x){
  res <- x - apply(x, 2, median)
  return(res)
}

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") %>%
  pull(cell_id)

data.plot <- data.plot[,tmp.select]
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)

#fviz_eig

eig <- factoextra::get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(beforeEPI.meta,by = "cell_id")

data.plot$Stage[1:3] <- "MII oocyte" 


Stage <- c("MII oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast")

data.plot$Stage <- factor(data.plot$Stage,levels = Stage)

###edit embryo id
data.plot$embryo_id[1:3] <- 1:3
data.plot$embryo_id[4:7] <- 1:4

data.plot <- data.plot %>%
  mutate_cond(Stage == "early2cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c("0r",1,2,3),to = 1:4)) %>%
  mutate_cond(Stage == "mid2cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c("0r",3:7),to = 1:6)) %>%
  mutate_cond(Stage == "late2cell",
              embryo_id = plyr::mapvalues(embryo_id,from = 5:9,to = 1:5)) %>%
  mutate_cond(Stage == "8cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c(1,2,5,8),to = c(1,2,3,4))) %>%
  mutate_cond(Stage == "16cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c(1,4,5,6),to = c(1,2,3,4))) %>%
  mutate_cond(Stage == "earlyblast",
              embryo_id = plyr::mapvalues(embryo_id,from = 2:4,to = c(1,2,3)))



ggplot(data.plot,aes(PC1,PC2,color=Stage))+
  geom_point(size=3)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(24,25,16,17,18,19))+
  theme_cowplot()

ggsave(filename = myFileName(prefix = "res/fig/QC_replicate_deng_fig1A_v1",suffix = ".jpg"),
       width = 8,height = 8,dpi = 300)

### more precise
ggplot(data.plot,aes(PC1,PC2,color=Stage,shape=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  theme_cowplot()

ggsave(filename = myFileName(prefix = "res/fig/QC_replicate_deng_fig1A_v2",suffix = ".jpg"),
       width = 8,height = 8,dpi = 300)

######Next will putinto the main figure
######-------------2.3.2  add xiewei's data-----------------

######-------------2.3.2.1 boxplot-----------

######pseudo-bulk 

Stage <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast")

data.plot.1 <- beforeEPI.exp.mean[,Stage]
tmp.names.1 <- colnames(data.plot.1)
tmp.names.2 <- rownames(data.plot.1)
data.plot.1 <- preprocessCore::normalize.quantiles(as.matrix(data.plot.1))
colnames(data.plot.1) <- tmp.names.1
rownames(data.plot.1) <- tmp.names.2

data.plot.1 <- data.plot.1 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "gene.exp",-Gene) 


data.plot.1$Stage <- factor(data.plot.1$Stage,Stage)

tmp.stage <-  c("MII oocyte","zygote","early 2-cell","mid 2-cell","late 2-cell",
                "4-cell","8-cell","16-cell","early blastocyst","mid blastocyst",  
                "late blastocyst")
p <- myBoxplot_advanced(data = data.plot.1,
                        x = "Stage",
                        y = "gene.exp",
                        fontsize = 30,
                        axislinesize = 2,
                        fill = "Stage",color = NULL,
                        mycolor = test.color.2(length(Stage)),
                        ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("gene.exp.mean")+
  xlab(NULL)+
  theme(legend.position = "none")
myggsave(p = p,prefix = "res/fig/figS1_QC_pseudobulk_boxplot",suffix = ".jpg",width = 8,height = 8,dpi = 350)


tmp.stage <- c("MII_oocyte","zygote","early2cell","late2cell",
               "4cell","8cell")

data.plot.2 <- xiewei.data.exp[,tmp.stage]
tmp.names.1 <- colnames(data.plot.2)
tmp.names.2 <- rownames(data.plot.2)
data.plot.2 <- preprocessCore::normalize.quantiles(as.matrix(data.plot.2))
colnames(data.plot.2) <- tmp.names.1
rownames(data.plot.2) <- tmp.names.2


data.plot.2 <- data.plot.2 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp,-Gene) %>%
  mutate(group = "bulk")


data.plot.2$Stage <- factor(data.plot.2$Stage,levels = tmp.stage)

tmp.stage <- c("MII oocyte","zygote","early 2-cell","late 2-cell",
               "4-cell","8-cell")
p <- myBoxplot_advanced(data = data.plot.2,
                        x = "Stage",
                        y = "gene.exp",
                        fill = "Stage",
                        color = NULL,fontsize = 30,axislinesize = 2,
                        mycolor = test.color.2(length(Stage)),
                        ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("log2(RPKM+1)")+
  xlab(NULL)+
  theme(legend.position = "none")
p
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot",suffix = ".jpg",width = 8,height = 8,dpi = 350)

#####-----2.3.2.2 PCA---------

gene.use <- intersect(rownames(xiewei.data.exp),rownames(beforeEPI.exp.mean))
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") %>%
  pull(cell_id)
data.plot <- data.plot[gene.use,tmp.select]

tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")
data.plot <- cbind(xiewei.data.exp[gene.use,tmp.select],data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####add meta data for xieweidata
tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")
ttttt <- data.frame(cell_id=tmp.select,
                    cell_type=tmp.select,
                    embryo_id="x",
                    embryo_cell_id = tmp.select,
                    Stage=tmp.select,
                    stringsAsFactors = F)

beforeEPI.meta$embryo_id[1:3] <- 1:3
beforeEPI.meta$embryo_id[4:7] <- 1:4

# beforeEPI.meta <- beforeEPI.meta %>%
#   mutate_cond(Stage == "early2cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c("0r",1,2,3),to = 1:4)) %>%
#   mutate_cond(Stage == "mid2cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c("0r",3:7),to = 1:6)) %>%
#   mutate_cond(Stage == "late2cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = 5:9,to = 1:5)) %>%
#   mutate_cond(Stage == "8cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c(1,2,5,8),to = c(1,2,3,4))) %>%
#   mutate_cond(Stage == "16cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c(1,4,5,6),to = c(1,2,3,4))) %>%
#   mutate_cond(Stage == "earlyblast",
#               embryo_id = plyr::mapvalues(embryo_id,from = 2:4,to = c(1,2,3)))

ttttt$data.type <- "RNA-seq"
beforeEPI.meta$data.type <- "scRNA-seq"
tmp.meta <- rbind(ttttt,beforeEPI.meta)


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")

data.plot$Stage[c(1,7:9)] <- "MII_oocyte" 

Stage <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast")

Stage.1 <- c("MII oocyte","zygote",
             "early 2-cell","mid 2-cell","late 2-cell",
             "4-cell","8-cell","16-cell",
             "early blastocyst","mid blastocyst","late blastocyst")
data.plot$Stage <- plyr::mapvalues(data.plot$Stage,from = Stage,to = Stage.1)

unique(data.plot$Stage)

data.plot$Stage <- factor(data.plot$Stage,levels = Stage.1)


p <- ggplot(data.plot,aes(PC1,PC2,color=Stage,shape=data.type))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 24)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 2 ),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
myggsave(p,prefix = "res/fig/Fig1_QC_pca",suffix = ".jpg",width = 8,height = 6,dpi = 350)


###-------2.4 merge t-SNE----------


seu <- readRDS(file = "res/R/early.scRNAseq.seurat_20210130.Rds")
# tmp.colors <- c("#FCB03B","#0F6937","#67BD45",
#   "#169192","#B09DBA","#954496",
#   "#4E1550","#8EAF3D","#F05153",
#   "#5358A5","#79A7AA","#62BDA6","#4098B6","#4173B3","#5E4FA2","#000000")

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst",
                  "E5.25","E5.5","E6.25","E6.5","fibroblast")
p <- TSNEPlot(seu,label=F,group="Stage",pt.size=1.5)+
  ggtitle(NULL)+
  scale_color_manual(values = col.spectral(length(tmp.levels.2)))+
  theme_cowplot(font_size = 24)+
  theme(axis.line = element_line(size = 2),aspect.ratio = 1.5)
p
ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


####--------3. Ligand, Receptor expression-----------
#####-------3.1 show Ligand, Receptor gene expression landscape--------
#####--------3.1.1 prepare for Ligand, Receptor database-----------
LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",stringsAsFactors = F)
LRpairs <- LRpairs.df$lr_pair
Lgenelist <- LRpairs.df$ligand_gene_symbol
Rgenelist <- LRpairs.df$receptor_gene_symbol 

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")

gene_symbols <- rownames(data.merge.exp)
colnames(data.merge.meta)[1] <- "Cell"
l.remove <- setdiff(Lgenelist,gene_symbols)
r.remove <- setdiff(Rgenelist,gene_symbols)
index.remove <- c(which(Lgenelist %in% l.remove),which(Rgenelist %in% r.remove))
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
saveRDS(LRpairs,file = myFileName(prefix = "res/R/cytotalkdb_LRpairs",suffix = ".rds"))

#####------3.1.2 prepare data for plot-------
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")
colnames(data.merge.meta)[1] <- "Cell"

tmp.row.names <- rownames(data.merge.exp)
tmp.col.names <- colnames(data.merge.exp)

data.merge.exp <- preprocessCore::normalize.quantiles(as.matrix(data.merge.exp))
rownames(data.merge.exp) <- tmp.row.names
colnames(data.merge.exp) <- tmp.col.names

data.merge.exp.L <- data.merge.exp[Lgene,]
data.merge.exp.R <- data.merge.exp[Rgene,]

data.plot.L <- data.merge.exp.L %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  gather(key = Cell,gene.exp,-Gene) %>%
  full_join(y=data.merge.meta,by = "Cell")
data.plot.L$group <- "Lgene"

#### this was used as test
data.plot.L %>%
  rowid_to_column() %>%
  dplyr::filter(is.na(gene.exp)) %>%
  head()

data.plot.R <- data.merge.exp.R %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  gather(key = Cell,gene.exp,-Gene) %>%
  full_join(y=data.merge.meta,by = "Cell")
data.plot.R$group <- "Rgene"

data.plot <-  data.merge.exp %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  gather(key = Cell,gene.exp,-Gene) %>%
  full_join(y=data.merge.meta,by = "Cell")
data.plot$group <- "Background"
nrow(data.plot)
data.plot <- rbind(data.plot.L,data.plot.R,data.plot)
nrow(data.plot)
cat(paste0("\"",unique(data.plot$Stage),"\""),sep = ",")

### prepare for data.plot
data.plot$Stage <- factor(data.plot$Stage,
                          levels = c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                                     "4cell","8cell","16cell","earlyblast","midblast","lateblast",
                                     "E5.25","E5.5","E6.25","E6.5","fibroblast"))

saveRDS(data.plot,file = myFileName(prefix = "res/R/data_merge_plot",suffix = ".rds"))



#####---------3.1.3 plot combine geom_boxplot----------------


#####---------3.1.3.1 load data and checkout----------
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_20211005.rds")
###checkout the group assignment
data.plot %>%
  group_by(Cell,group) %>%
  summarise(n=n()) %>%
  head()
###check normalization
summary(data.plot$gene.exp)
###check na
data.plot %>%
  rowid_to_column() %>%
  dplyr::filter(is.na(gene.exp)) %>%
  head()
###
levels(data.plot$Stage)
head(data.plot)
data.plot$group <- factor(data.plot$group,levels = c("Lgene","Rgene","Background"))


#####-------3.1.3.2 raw gene.exp boxplot-------

tmp.data.plot <- data.plot
mypalette <- c("#ffc000","#00b050","darkgrey")
head(tmp.data.plot)
scales::show_col(mypalette)
names(mypalette) <- c("Ligand","Receptor","Background")
x <- "Stage"
y <- "gene.exp"
head(tmp.data.plot)
p1 <- ggplot()+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp=gene.exp),
               aes(Stage,gene.exp,fill=names(mypalette)[3],color=names(mypalette)[3]),
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp),
               aes(Stage,gene.exp,fill=names(mypalette)[1],color=names(mypalette)[1]),
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp),
               aes(Stage,gene.exp),
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp),
               aes(Stage,gene.exp,group=""),
               geom = "line", 
               color="black",
               fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp=-gene.exp),
               aes(Stage,gene.exp,fill=names(mypalette)[3],color=names(mypalette)[3]),
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp),
               aes(Stage,-gene.exp,fill=names(mypalette)[1],color=names(mypalette)[1]),
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp),
               aes(Stage,-gene.exp),
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp),
               aes(Stage,-gene.exp,group=""),
               geom = "line", 
               color="black",
               fun = mean)+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(limits = c(-6,6),
                     breaks = seq(-6,6,2),
                     labels = as.character(c(-seq(-6,-2,by = 2),0,seq(2,6,by = 2))))+
  ylab("log2(RPKM+1)")+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))

myggsave(p = p1,prefix = "res/fig/Fig2_landscape_gene_exp",suffix = ".jpg",width = 8,height = 6,dpi=350)


####--------------3.2 gene mean expression changes------------------

####-------------3.2.1 calc mean expression-------------
####mean expression
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_20211005.rds")

#####Background
data.plot.mean.Background <- data.plot %>%
  dplyr::filter(group=="Background") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()
#####Ligand
data.plot.mean.Ligand <- data.plot %>%
  dplyr::filter(group=="Lgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()

#####Receptor
data.plot.mean.Receptor <- data.plot %>%
  dplyr::filter(group=="Rgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()



data.plot.mean.Background$group <- "Background"
data.plot.mean.Ligand$group <- "Lgene"
data.plot.mean.Receptor$group <- "Rgene"

data.plot.mean <- rbind(data.plot.mean.Background,
                        data.plot.mean.Ligand,
                        data.plot.mean.Receptor)

data.plot.mean$group <- factor(data.plot.mean$group,
                               levels = c("Lgene","Rgene","Background"))


#####------3.2.2 plot-------


#####-----3.2.2.1 mean--------
tmp.data.plot <- data.plot.mean
levels(tmp.data.plot$group)
mypalette <- c("#ffc000","#00b050","lightgrey")

p2 <- ggplot(data = tmp.data.plot,aes(Stage,gene.exp.mean))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.mean=gene.exp.mean),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean),
               fill = mypalette[1], 
               color=mypalette[1],
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), 
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), 
               geom = "line", 
               color="black",
               aes(group=""),
               fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.mean=-gene.exp.mean),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.mean=-gene.exp.mean),
               fill = mypalette[2], color=mypalette[2],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.mean=-gene.exp.mean), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.mean=-gene.exp.mean), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_hline(yintercept = 0,linetype="solid")+
  ylab("gene.exp.mean")+
  scale_y_continuous(limits = c(-6,6),
                     breaks = seq(-6,6,2),
                     labels = as.character(c(-seq(-6,-2,by = 2),0,seq(2,6,by = 2))))+
  scale_x_discrete(labels=tmp.stage <- c("MII oocyte","zygote",
                                         "early 2-cell","mid 2-cell","late 2-cell",
                                         "4-cell","8-cell","16-cell",
                                         "early blastocyst","mid blastocyst","late blastocyst",
                                         "E5.25","E5.5","E6.25","E6.5","fibroblast"))+
  theme_cowplot(font_size  = 28)+
  theme(legend.position = "right",
        axis.ticks.y = element_line(size = 2),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))


p2 <- p2 + ylab("gene.exp.mean")+xlab(NULL)
p2
myggsave(p = p2,prefix = "res/fig/Fig2_landscape_gene_mean_exp",
         suffix = ".jpg",
         width = 8,height = 6,dpi=350)
myggsave(p = p2,prefix = "res/fig/Fig2_landscape_gene_mean_exp",
         suffix = ".pdf",
         width = 8,height = 6)

p2+coord_polar()+theme_map()

p3 <- ggplot(data = tmp.data.plot,aes(Stage,gene.exp.sd))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.sd=gene.exp.sd),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.sd=gene.exp.sd),
               fill = mypalette[1], color=mypalette[1],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.sd=gene.exp.sd), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.sd=gene.exp.sd), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.sd=-gene.exp.sd),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.sd=-gene.exp.sd),
               fill = mypalette[2], color=mypalette[2],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.sd=-gene.exp.sd), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.sd=-gene.exp.sd), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(limits = c(-5,5),
                     breaks = seq(-5,5,1),
                     labels = as.character(c(-seq(-5,-1,by = 1),0,seq(1,5,by = 1))))+
  ylab("gene.exp.sd")+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))

myggsave(p = p3,prefix = "res/fig/FigS3_landscape_gene_sd",suffix = ".jpg",width = 8,height = 6,dpi=350)




p4 <- ggplot(data = tmp.data.plot,aes(Stage,gene.exp.cv))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.cv=gene.exp.cv),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.cv=gene.exp.cv),
               fill = mypalette[1], color=mypalette[1],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.cv=gene.exp.cv), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.cv=gene.exp.cv), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.cv=-gene.exp.cv),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.cv=-gene.exp.cv),
               fill = mypalette[2], color=mypalette[2],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.cv=-gene.exp.cv), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.cv=-gene.exp.cv), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(limits = c(-30,30),
                     breaks = seq(-30,30,10),
                     labels = as.character(c(-seq(-30,-10,by = 10),0,seq(10,30,by = 10))))+
  ylab("gene.exp.cv")+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
#p3 
myggsave(p = p4,prefix = "res/fig/FigS3_gene_cv",suffix = ".jpg",width = 8,height = 6,dpi=300)

#####---------3.2.2.2 rank---------------

LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

data.plot.rank.Background <- data.plot.mean.Background %>%
  group_by(Stage) %>%
  mutate(rank=row_number(gene.exp.mean)) %>%
  ungroup() %>%
  dplyr::select(-group,-gene.exp.mean,-gene.exp.sd,-gene.exp.cv) %>%
  spread(key = "Stage",value = "rank") %>%
  column_to_rownames("Gene")

data.plot.rank.L <- data.plot.rank.Background[Lgene,] %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "rank",-Gene) %>%
  mutate(group="Lgene")

data.plot.rank.R <- data.plot.rank.Background[Rgene,] %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "rank",-Gene) %>%
  mutate(group="Rgene")

data.plot.rank <- data.plot.rank.Background %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "rank",-Gene) %>%
  mutate(group="Background")

data.plot.rank <- rbind(data.plot.rank,data.plot.rank.L,data.plot.rank.R)

data.plot.rank$Stage <- factor(data.plot.rank$Stage,
                               levels = c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                                          "4cell","8cell","16cell","earlyblast","midblast","lateblast",
                                          "E5.25","E5.5","E6.25","E6.5","fibroblast"))

data.plot.rank$group <- factor(data.plot.rank$group,
                               levels = c("Lgene","Rgene","Background"))



tmp.data.plot <- data.plot.rank
head(tmp.data.plot)
mypalette <- c("#ffc000","#00b050","lightgrey")



tmp.data.plot.bar <- subset(tmp.data.plot,group == "Background") %>%
  mutate(rank=rank) %>%
  group_by(Stage) %>%
  summarise(rank=n()) %>%
  ungroup()

levels(tmp.data.plot$Stage)[1:11]

tmp.stage <- c("MII oocyte","zygote",
  "early 2-cell","mid 2-cell","late 2-cell",
  "4-cell","8-cell","16-cell",
  "early blastocyst","mid blastocyst","late blastocyst",
  "E5.25","E5.5","E6.25","E6.5","fibroblast")
p5 <- ggplot(data = tmp.data.plot,aes(Stage,rank))+
  # geom_bar(data = tmp.data.plot.bar,
  #          mapping = aes(x = Stage,y = rank),fill = mypalette[3], color=mypalette[3],
  #          stat = "identity",width = 0.8)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(rank=rank), geom = "point",shape = 15, fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(rank=rank), geom = "line", color=mypalette[1],aes(group=""),fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(rank=rank), geom = "point",shape = 15, fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(rank=rank), geom = "line", color=mypalette[2],aes(group=""),fun = mean)+
  geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
  scale_y_continuous(breaks = seq(6500,8500,500))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("rank")+
  xlab(NULL)+
  theme_cowplot(font_size  = 28)+
  theme(legend.position = "right",
        axis.ticks = element_line(size = 2),
        axis.line.x = element_line( size = 2 ),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 35))
p5
myggsave(p = p5,prefix = "res/fig/Fig2_gene_exp.mean_rank",suffix = ".jpg",width = 8,height = 6,dpi=350)


####
x <- tmp.data.plot %>%
  dplyr::filter(group=="Lgene") %>%
  group_by(Stage) %>%
  summarise(ttt=mean(rank)) %>%
  pull(ttt)
y <- tmp.data.plot %>%
  dplyr::filter(group=="Rgene") %>%
  group_by(Stage) %>%
  summarise(ttt=mean(rank)) %>%
  pull(ttt)
tmp.res <- cor.test(x[1:11],y[1:11])
tmp.res$p.value
tmp.stage


p <- p5+annotate("text",
                 x = "4cell",
                 y = 8500,
                 label=paste0("PCC=",round(cor(x[1:11],y[1:11]),4),"\n",
                              "p=",round(cor.test(x[1:11],y[1:11])$p.value,4)),
                 size=10)
p
myggsave(p = p,prefix = "res/fig/Fig2_gene_exp.mean_rank",suffix = ".jpg",
         width = 8,height = 6,dpi=350)



##-------8.4 visualize gene expression---------

###load data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)
Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","fibroblast")

LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20201214.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast","lateblast")

tmp.cell_id <- beforeEPI.meta %>%
  dplyr::filter(Stage %in% Stage.level) %>%
  pull(cell_id)

### corplot
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
hist(tmp.res.cor.1$PCC)
cor.test(rnorm(100),rnorm(100))
test.color <- colorRampPalette(RColorBrewer::brewer.pal(name = "OrRd",n=9))
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
Stage.level <-  c("Oocyte","zygote","early2cell",
                  "mid2cell","late2cell","4cell",
                  "8cell","16cell","earlyblast",
                  "midblast","lateblast")
paper.color <- c("#FCB03B","#0F6937","#67BD45",
                 "#169192","#B09DBA","#954496",
                 "#4E1550","#8EAF3D","#F05153",
                 "#5358A5","#79A7AA")
scales::show_col(paper.color)
myCorplot <- function(mat=beforeEPI.exp[,tmp.cell_id],
                      mat.meta=beforeEPI.meta,
                      tmp.Lgene=tmp.res.cor.1$Lgene[1],
                      tmp.Rgene=tmp.res.cor.1$Rgene[1],
                      PCC=tmp.res.cor.1$PCC[1],
                      Stage.level=Stage.level,
                      color=test.color(length(Stage.level))){
  
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(t(mat.L),
                       t(mat.R),
                       stringsAsFactors = F) %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    mutate(Stage=factor(Stage,levels = Stage.level))
  
  p<- ggplot(data = res.df,
             aes_string(x=tmp.Lgene,y=tmp.Rgene,color="Stage"))+
    geom_point(size=8)+
    scale_color_manual(values = color)+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 15)+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 15 ) ,
          axis.text = element_text(size = 15),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ))
  return(p)
}

tmp.idx <- which(tmp.res.cor.1$Lgene %in% "Fgf4")

plot.list <- list()
for(ii in tmp.idx){
  p <- myCorplot(mat=beforeEPI.exp[,tmp.cell_id],
                 mat.meta=beforeEPI.meta,
                 tmp.Lgene=tmp.res.cor.1$Lgene[ii],
                 tmp.Rgene=tmp.res.cor.1$Rgene[ii],
                 PCC=tmp.res.cor.1$PCC[ii],
                 Stage.level=Stage.level,
                 color=test.color.2(length(Stage.level)))
  plot.list <- c(plot.list,list(p))
}
plot_grid(plotlist = plot.list)





####scatter plot

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")

Stage.level <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast",  
                  "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")
tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")

tmp.cell_id <- data.merge.meta %>%
  dplyr::filter(Stage %in% Stage.level) %>%
  pull(cell_id)


myScatterPlot <- function(mat = data.merge.exp[,tmp.cell_id],
                          mat.meta = data.merge.meta,
                          tmp.Lgene = tmp.data.plot$Lgene[nrow(tmp.data.plot)],
                          tmp.Rgene = tmp.data.plot$Rgene[nrow(tmp.data.plot)],
                          PCC = tmp.data.plot$PCC[nrow(tmp.data.plot)]){
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(t(mat.L),
                       t(mat.R),
                       stringsAsFactors = F) %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    mutate(Stage=factor(Stage,levels = Stage.level)) %>%
    gather(key = "gene",value = "gene.exp",-cell_id,-Stage) %>%
    mutate(gene=factor(gene,levels = c(tmp.Lgene,tmp.Rgene)))
  
  
  
  p <- ggplot(data = res.df,aes_string(x="Stage",y="gene.exp",colour="gene"))+
    geom_point(alpha=0.1,size=1)+
    geom_smooth(method = "gam",aes_string(group="gene"),se = F)+
    #scale_color_manual(values = c("#E41A1C","#377EB8"))+
    scale_color_manual(values = c("#ffc000","#00b050"))+
    ylab("gene.exp")+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 28)+
    theme(legend.position = "top",
          legend.justification = "center",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 28 ) ,
          axis.text = element_text(size = 28),
          axis.line.x = element_line( size = 2 ),
          axis.line.y = element_line( size = 2 ),
          axis.text.x = element_text(hjust = 0.9,angle = 45))
  return(p)
}
#####structure!!
# p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
#                    mat.meta = data.merge.meta,
#                    tmp.Lgene = "Ctcf",
#                    tmp.Rgene = "Gtf3c2",
#                    PCC = "Ctcf")+
#   ggtitle("Ctcf")+
#   scale_x_discrete(label=tmp.stage)
# 
# p

####"Bmp6_Bmpr1b"
tmp.idx <- 4
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = sprintf("%0.3f",tmp.res.cor.1$PCC[tmp.idx]))+
  scale_x_discrete(label=tmp.stage)+
  geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
  theme(axis.ticks = element_line(size = 2))+
  xlab(NULL)
p

ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

ggsave(p,filename = myFileName(prefix = paste0("res/fig/fig2_eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".pdf"),
       width = 8,height = 8)

tmp.idx <- 1072
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = sprintf("%0.3f",tmp.res.cor.1$PCC[tmp.idx]))+
  scale_x_discrete(label=tmp.stage)+
  geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
  theme(axis.ticks = element_line(size = 2))+
  xlab(NULL)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

ggsave(p,filename = myFileName(prefix = paste0("res/fig/fig2_eLR_",
                                               tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".pdf"),
       width = 8,height = 8)


#####Igf2-Igf2r

tmp.idx <- 983
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = sprintf("%0.3f",tmp.res.cor.1$PCC[tmp.idx]))+
  scale_x_discrete(label=tmp.stage)+
  geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
  theme(axis.ticks = element_line(size = 2))+
  xlab(NULL)
p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/fig2_eLR_",
                                               tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".pdf"),
       width = 8,height = 8)


tmp.idx <- 978
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = sprintf("%0.3f",tmp.res.cor.1$PCC[tmp.idx]))+
  scale_x_discrete(label=tmp.stage)+
  geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
  theme(axis.ticks = element_line(size = 2))+
  xlab(NULL)
p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/fig2_eLR_",
                                               tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".pdf"),
       width = 8,height = 8)



tmp.idx <- 264
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = sprintf("%0.3f",tmp.res.cor.1$PCC[tmp.idx]))+
  scale_x_discrete(label=tmp.stage)+
  geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
  theme(axis.ticks = element_line(size = 2))+
  xlab(NULL)
p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/fig2_eLR_",
                                               tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".pdf"),
       width = 8,height = 8)

####
tmp.idx <- 1007
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = sprintf("%0.3f",tmp.res.cor.1$PCC[tmp.idx]))+
  scale_x_discrete(label=tmp.stage)+
  geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
  theme(axis.ticks = element_line(size = 2))+
  xlab(NULL)
p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/fig2_eLR_",
                                               tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".pdf"),
       width = 8,height = 8)




