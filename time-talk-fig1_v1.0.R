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
library(rtracklayer)
library(GenomicFeatures)
library(patchwork)
source(file = "code/myUtils.R")

# ###-----keep you on R------
# for(ii in 1:1000000){
#   cat(paste("round",ii),sep = "\n")
#   Sys.sleep(ii/10)
# }



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

beforeEPI.exp.mean$ICM <- beforeEPI.exp.mean %>%
  dplyr::select(c("earlyblast","midblast","lateblast")) %>%
  rowMeans()


####-------revised version-----------
Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
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
         width = 16,height = 16,dpi=350)


#####-----pheatmap-----

Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
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

tmp.stage  <- c("MII oocyte","zygote","early 2-cell","late 2-cell","4-cell","8-cell","ICM")

png(filename = myFileName(prefix = "res/fig/Fig1_QC_pheatmap",suffix = ".png"),
     width = 8,height = 8,res = 350,units = "in")
exprTable_t <- as.data.frame(t(res.mat))
col_dist = dist(exprTable_t)
hclust_1 <- hclust(col_dist)
manual_order = rownames(res.mat)
#dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))))
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))), agglo.FUN = max)
col_cluster <- as.hclust(dend)

ComplexHeatmap::pheatmap(res.mat,fontsize = 20,
         show_rownames = T,
         display_numbers = T,
         number_color = "white",
         cluster_rows = col_cluster,
         cluster_cols = col_cluster,
         fontsize_number = 22,
         name = "PCC",
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

######------2.3.2.1 pseudo-bulk---------

######normalized
# Stage <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
#            "4cell","8cell","16cell","earlyblast","midblast",
#            "lateblast")

Stage <-  c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
             "4cell","8cell","16cell","earlyblast","midblast",
             "lateblast","E5.25","E5.5","E6.25","E6.5")

seu <- readRDS("res/R/early_gastrulation_20210115.rds")
seu@assays$RNA@counts <- 2^(seu@assays$RNA@data)-1
Idents(seu) <- "Stage"
tmp.post.implantation <- AverageExpression(seu,slot = "counts")$RNA
tmp.post.implantation <- log2(tmp.post.implantation+1)

gene_symbols <- intersect(rownames(tmp.post.implantation),
                          rownames(beforeEPI.exp.mean))
sc.pseudobulk <- cbind(beforeEPI.exp.mean[gene_symbols,],
                       tmp.post.implantation[gene_symbols,])


#####---------quantile normalized---------

data.plot.1 <- sc.pseudobulk[,Stage]
colnames(sc.pseudobulk)
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
                "late blastocyst","E5.25","E5.5","E6.25","E6.5")
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
  theme(legend.position = "none",axis.ticks = element_line(size = 2))
myggsave(p,
         prefix = "res/fig/figS1_QC_pseudobulk_boxplot",
         suffix = ".png",
         width = 12,height = 9,dpi = 350)


#####------raw unquantile normlized------

data.plot.1 <- sc.pseudobulk[,Stage]
# colnames(sc.pseudobulk)
# tmp.names.1 <- colnames(data.plot.1)
# tmp.names.2 <- rownames(data.plot.1)
# data.plot.1 <- preprocessCore::normalize.quantiles(as.matrix(data.plot.1))
# colnames(data.plot.1) <- tmp.names.1
# rownames(data.plot.1) <- tmp.names.2

data.plot.1 <- data.plot.1 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "gene.exp",-Gene) 


data.plot.1$Stage <- factor(data.plot.1$Stage,Stage)

tmp.stage <-  c("MII oocyte","zygote","early 2-cell","mid 2-cell","late 2-cell",
                "4-cell","8-cell","16-cell","early blastocyst","mid blastocyst",  
                "late blastocyst","E5.25","E5.5","E6.25","E6.5")
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
  theme(legend.position = "none",axis.ticks = element_line(size = 2))
p
myggsave(p = p,
         prefix = "res/fig/figS1_QC_pseudobulk_boxplot_raw",
         suffix = ".png",
         width = 12,height = 9,dpi = 350)




####----2.3.3.2 xieweidata -----

#####unormlized
tmp.stage <- c("MII_oocyte","zygote","early2cell","late2cell",
               "4cell","8cell","ICM")

data.plot.2 <- xiewei.data.exp[,tmp.stage]
# tmp.names.1 <- colnames(data.plot.2)
# tmp.names.2 <- rownames(data.plot.2)
# data.plot.2 <- preprocessCore::normalize.quantiles(as.matrix(data.plot.2))
# colnames(data.plot.2) <- tmp.names.1
# rownames(data.plot.2) <- tmp.names.2


data.plot.2 <- data.plot.2 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp,-Gene) %>%
  mutate(group = "bulk")


data.plot.2$Stage <- factor(data.plot.2$Stage,levels = tmp.stage)

tmp.stage <- c("MII oocyte","zygote","early 2-cell","late 2-cell",
               "4-cell","8-cell","ICM")
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
  theme(legend.position = "none",axis.ticks = element_line(size = 2))
p
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot",
         suffix = ".jpg",width = 8,height = 8,dpi = 350)





####normalized data

tmp.stage <- c("MII_oocyte","zygote","early2cell","late2cell",
               "4cell","8cell","ICM")

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
               "4-cell","8-cell","ICM")
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
  theme(legend.position = "none",axis.ticks = element_line(size = 2))
p
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot_normalized",
         suffix = ".jpg",width = 8,height = 8,dpi = 350)

#####-----2.3.2.2 PCA---------

gene.use <- intersect(rownames(xiewei.data.exp),rownames(beforeEPI.exp.mean))
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") %>%
  pull(cell_id)
data.plot <- data.plot[gene.use,tmp.select]

tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
data.plot <- cbind(xiewei.data.exp[gene.use,tmp.select],data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####add meta data for xieweidata
tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
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

data.plot$Stage[c(1,8:10)] <- "MII_oocyte" 

Stage <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","ICM","earlyblast","midblast",  
           "lateblast")

Stage.1 <- c("MII oocyte","zygote","early 2-cell",
             "mid 2-cell","late 2-cell","4-cell",
             "8-cell","16-cell","ICM",
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
                               "#4E1550","#8EAF3D","grey",
                               "#F05153","#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 24)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 2 ),
        axis.line.y = element_line( size = 2 ),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
p
myggsave(p,prefix = "res/fig/Fig1_QC_pca",suffix = ".jpg",width = 8,height = 6,dpi = 350)
myggsave(p,prefix = "res/fig/Fig1_QC_pca",suffix = ".pdf",width = 8,height = 6)

###-------2.4 merge t-SNE----------


seu <- readRDS(file = "res/R/early.scRNAseq.seurat_20211218.rds")
# tmp.colors <- c("#FCB03B","#0F6937","#67BD45",
#   "#169192","#B09DBA","#954496",
#   "#4E1550","#8EAF3D","#F05153",
#   "#5358A5","#79A7AA","#62BDA6","#4098B6","#4173B3","#5E4FA2","#000000")
seu <- subset(seu,Stage!="fibroblast")
colnames(seu[[]])

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst",
                  "E5.25","E5.5","E6.25","E6.5")
p <- TSNEPlot(seu,label=F,group="Stage",pt.size=1.5)+
  ggtitle(NULL)+
  scale_color_manual(values = col.spectral(length(tmp.levels.2)))+
  theme_cowplot()+
  geom_path(data = data.frame(x=c(-45,-35),y=c(-35,-35)),
            size = 1,
            mapping = aes(x,y),linejoin = "bevel", lineend = "round",
            arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches")))+
  geom_path(data = data.frame(x=c(-45,-45),y=c(-35,-25)),
            size = 1,
            mapping = aes(x,y),linejoin = "bevel", lineend = "round",
            arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches")))+
  annotate("text",x = -40,y = -40,label="tSNE1",fontface="bold",size=8)+
  annotate("text",x = -50,y = -30,label="tSNE2",fontface="bold",size=8,angle=90)+
  NoAxes()+
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 24))
p
ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne",suffix = ".pdf"),
       width = 8,height = 8)

p1 <- FeaturePlot(seu,features = c("Sox2","Pou5f1","Nanog","Gata6","Zscan4c","Bmp4")) & 
  scale_color_gradientn(colors = rdwhbu(100)) &
  theme(axis.line = element_line(size = 2),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 2))
p1
ggsave(filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".png"),
       width = 8,height = 9,dpi = 350)


##### show lineage

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst",
                  "Epiblast(EPI),","Ectoderm(ExE)",
                  "Visceral Endoderm (VE)","fibroblast")
p <- TSNEPlot(seu,label=F,group="cell_type",pt.size=1.5)+
  ggtitle(NULL)+
  scale_color_manual(values = col.spectral(length(tmp.levels.2)))+
  theme_cowplot()+
  geom_path(data = data.frame(x=c(-45,-35),y=c(-35,-35)),
            size = 1,
            mapping = aes(x,y),linejoin = "bevel", 
            lineend = "round",
            arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches")))+
  geom_path(data = data.frame(x=c(-45,-45),y=c(-35,-25)),
            size = 1,
            mapping = aes(x,y),linejoin = "bevel", lineend = "round",
            arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches")))+
  annotate("text",x = -40,y = -40,label="tSNE1",
           fontface="bold",size=6)+
  annotate("text",x = -50,y = -30,label="tSNE2",
           fontface="bold",size=6,angle=90)+
  NoAxes()+
  theme(aspect.ratio = 4/3,
        legend.text = element_text(size = 24))
p
ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne_by_lineage",
                             suffix = ".jpg"),
       width = 8,height = 6,dpi = 350)

ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne_by_lineage",
                             suffix = ".pdf"),
       width = 8,height = 8)

p1 <- FeaturePlot(seu,features = c("Sox2","Pou5f1",
                                   "Nanog","Gata6",
                                   "Zscan4c","Bmp4","Amn")) +
  plot_layout(ncol = 2) &
  scale_color_gradientn(colors = rdwhbu(100)) &
  theme(axis.line = element_line(size = 2),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 2))
?FeaturePlot
p1
ggsave(filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".png"),
       width = 8,height = 12,dpi = 350)

# p <- TSNEPlot(seu,label=F,group="Stage",pt.size=1.5)+
#   ggtitle(NULL)+
#   scale_color_manual(values = col.spectral(length(tmp.levels.2)))+
#   theme_cowplot(font_size = 24)+
#   theme(axis.line = element_line(size = 2),aspect.ratio = 1.5)+
#   NoLegend()+
#   NoAxes()
# p
# ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne_formodel",
#                              suffix = ".pdf"),
#        width = 8,height = 8,dpi = 350)


####------2.5 Run CCA-------
####beforeEPI
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
rownames(beforeEPI.meta) <- beforeEPI.meta$cell_id

seu.pre <- myCreateSeurat(data = beforeEPI.exp,
               meta.data = beforeEPI.meta)
seu.pre <- FindVariableFeatures(seu.pre,
                                selection.method = "vst",
                                nfeatures = 2000)

####Idents
seu.post <- readRDS("res/R/early_gastrulation_20210115.rds")
VariableFeaturePlot(seu.post)



###Integrate
seu.anchor <- FindIntegrationAnchors(object.list = list(seu.pre,seu.post),dims = 1:20)

seu.integrated <- IntegrateData(anchorset = seu.anchor, 
                                dims = 1:20)

dim(seu.integrated[["RNA"]]@counts)
DefaultAssay(seu.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
testAB.integrated <- seu.integrated
testAB.integrated <- ScaleData(testAB.integrated, 
                               features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, 
                            npcs = 50, 
                            verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, 
                                   dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, 
                                  resolution = 0.5)
testAB.integrated <- RunUMAP(testAB.integrated, 
                             dims = 1:30)
testAB.integrated <- RunTSNE(testAB.integrated, 
                             dims = 1:30)
TSNEPlot(testAB.integrated,group.by="Stage")
colnames(testAB.integrated[[]])



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

x <- "Stage"
y <- "gene.exp"

p1 <- ggplot(data = tmp.data.plot ,aes(Stage,gene.exp))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp=gene.exp),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp),
               fill = mypalette[1], 
               color=mypalette[1],
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp), 
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp), 
               geom = "line", 
               color="black",
               aes(group=""),
               fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp=-gene.exp),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp),
               fill = mypalette[2], 
               color=mypalette[2],
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp), 
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp), 
               geom = "line", 
               color="black",
               aes(group=""),
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
               fill = mypalette[1], color=mypalette[1],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), geom = "line", color="black",aes(group=""),fun = mean)+
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
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("gene.exp.mean")+
  scale_y_continuous(limits = c(-6,6),
                     breaks = seq(-6,6,2),
                     labels = as.character(c(-seq(-6,-2,by = 2),0,seq(2,6,by = 2))))+
  theme_cowplot(font_size  = 28)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 2 ),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))


p2 <- p2 + ylab("gene.exp.mean")+xlab(NULL)
myggsave(p = p2,prefix = "res/fig/Fig2_landscape_gene_mean_exp",
         suffix = ".jpg",
         width = 8,height = 6,dpi=350)




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
  #geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(breaks = seq(6500,8500,500))+
  ylab("rank")+
  xlab(NULL)+
  theme_cowplot(font_size  = 28)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 2 ),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 35))
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
cor(x[1:11],y[1:11])



p <- p5+annotate("text",x = "16cell",y = 8500,label=paste0("cor=",round(cor(x,y),4)),size=10)
myggsave(p = p,prefix = "res/fig/Fig2_gene_exp.mean_rank",suffix = ".jpg",width = 8,height = 6,dpi=350)



##-------4.visualize gene expression---------

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
  xlab(NULL)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),
                               suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 1072
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = sprintf("%0.3f",tmp.res.cor.1$PCC[tmp.idx]))+
  scale_x_discrete(label=tmp.stage)+
  xlab(NULL)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)



####-----------5. show track----------

####----------5.1 prepare genome--------

tmp.chrom <- read.delim(file = "database/mm9.chrom.sizes",header = F,stringsAsFactors = F)
head(tmp.chrom)
colnames(tmp.chrom) <- c("chrom","length")
tmp.metadata <- data.frame(name="Genome",value="mm9",stringsAsFactors = F)
####load txdb
mm9KG_txdb <- makeTxDbFromGFF(file = "D://Ubuntu/wlt/igenomes/Mus_musculus/UCSC/mm9/Annotation/genes.gtf",
                              chrominfo = tmp.chrom,
                              metadata = tmp.metadata)

#### the following function adapted from signac
txbd <- mm9KG_txdb
tx <- myGetGRangesFromsTxDb(txdb = mm9KG_txdb,standard.chromosomes = T,verbose = T)

tx[tx$gene_id=="Pou5f1",]

###--------5.2 test plot genes---------

tmp.dir <- "res/chipseq-process/bigwig/"
tmp.files <- list.files(path = tmp.dir,pattern = "bw",full.names = T)

tmp.levels <- c("MII_oocyte","PN5_zygote","2cell_early","2cell_late","4cell","8cell","ICM")
tmp.files <- paste0("res/chipseq-process/bigwig/GSE71434_",tmp.levels,"_K4me3.bed.bw")
seq_along(tmp.files)
###tmp.colors <- shaman.color(length(tmp.files))
tmp.colors <- c("#b4d1db","#6c1c1f","#4d311a","#d29156","#7f8a65","#6d82a3","#6e648a")
scales::show_col(tmp.colors)

####plot Zscan4c

tmp.region <- "chr7:11588091-11598898"

p_gene <- myGenePlot(annotation = tx,
                     region = tmp.region,
                     arrow_sbreaks = 600,
                     font_size = 18,label_size = 8)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.y = element_text(angle = 0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
p_gene



p_track_list <- lapply(seq_along(tmp.files),function(ii){
  p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                   bigwig = tmp.files[ii],
                   smooth = 100,
                   lognorm = T,
                   type = "coverage",
                   y_label = tmp.levels[ii],
                   fontsize=18,
                   track.color=tmp.colors[ii],
                   tmp.ylimits=c(0,0.01),
                   max.downsample = 3000,
                   downsample.rate = 0.1,
                   tmp.seed=42)+
    theme(axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1))
  p
  return(p)
})

p <- patchwork::wrap_plots(p_track_list,ncol = 1) +
  plot_annotation(title = tmp.region) & 
  theme(plot.title = element_text(hjust = 0.5,size = 18),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.y = element_text(size = 18,angle = 0,vjust = 0.5,hjust = -1),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.x = element_blank()) 
p1 <- p / p_gene

p_track_list <- lapply(seq_along(tmp.files),function(ii){
  p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                   bigwig = tmp.files[ii],
                   smooth = 100,
                   lognorm = T,
                   type = "coverage",
                   y_label = tmp.levels[ii],
                   fontsize=18,
                   track.color=tmp.colors[ii],
                   tmp.ylimits=c(0,5),
                   max.downsample = 3000,
                   downsample.rate = 0.1,
                   tmp.seed=42)+
    theme(axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1))
  return(p)
})
p_gene <- myGenePlot(annotation = tx,
                     region = tmp.region,
                     arrow_sbreaks = 4000,
                     font_size = 18,label_size = 8)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.y = element_text(angle = 0,vjust = 0.5,hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))

p <- patchwork::wrap_plots(p_track_list,ncol = 1) + 
  plot_annotation(title = tmp.region) & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        plot.title = element_text(size = 18,hjust = 0.5),
        axis.title.y = element_text(size = 18,angle = 0,vjust = 0.5,hjust = -1),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.x = element_blank())

p2 <- p / p_gene
p1 | (p2 & theme(axis.ticks.y = element_blank(),
                 axis.text.y  = element_blank(),
                 axis.title.y = element_blank()))

p1 | p2 

ggsave(filename = myFileName(prefix = "res/fig/fig5_epi_track_H3k4me3",suffix = ".png"),
       width = 16,height = 8,dpi = 350)

ggsave(p1,filename = myFileName(prefix = "res/fig/fig5_epi_track_fgf4_H3k4me3",suffix = ".png"),
       width = 16,height = 8,dpi = 350)

####----------5.3 write bedgraph to bigwig-------------

####test if work
# tmp.input <- "res/atac-seq-process/data/bedgraph/GSE66581_atac_ICM.bedGraph"
# tmp.output <- "res/tmp/test.bw"
# myBegraphToBigwig(tmp.input = tmp.input,
#                   tmp.output = tmp.output,
#                   region = "chr7:152042291-152056148",
#                   tmp.seqlength = seqlengths(tx))
### this is a bit time consuming

###---------5.3.1 ATAC-seq----------------------

tmp.files <- list.files(path = "res/chipseq-process/bedgraph/",full.names = T)
plan(multiprocess, workers = 8)
future_lapply(tmp.files,FUN = function(x){
  tmp.samples <- gsub(pattern = "res/atac-seq-process/data/bedgraph/GSE66581_atac_",replacement = "",x)
  tmp.samples <- gsub(pattern = ".bedGraph",replacement = ".bw",tmp.samples)
  tmp.output <- paste0("res/atac-seq-process/data/bigwig/GSE66581_atac_",tmp.samples)
  myBegraphToBigwig(tmp.input = x,
                    tmp.output = tmp.output,
                    region = NULL,
                    tmp.levels = seqlevels(tx),
                    tmp.seqlength = seqlengths(tx))
})
plan(sequential)
?plan
###----------5.3.2 H3K4me3 ----------------------
tmp.files <- list.files(path = "res/chipseq-process/bedgraph/",full.names = T)
plan(multiprocess, workers = length(tmp.files))
future_lapply(tmp.files,FUN = function(x){
  tmp.samples <- gsub(pattern = "res/chipseq-process/bedgraph/GSE71434_",replacement = "",x)
  tmp.samples <- gsub(pattern = ".bed",replacement = ".bw",tmp.samples)
  tmp.output <- paste0("res/chipseq-process/data_to_plot/H3K4me3/GSE71434_",tmp.samples)
  myBegraphToBigwig(tmp.input = x,
                    tmp.output = tmp.output,
                    region = NULL,
                    tmp.levels = seqlevels(tx),
                    tmp.seqlength = seqlengths(tx))
})
plan(sequential)
# region_data <- rtracklayer::import.bedGraph(
#   con = tmp.files[1],
#   which = as(region,"GRanges")
# )
# seqlevels(region_data) <- seqlevels(tx)
# seqlengths(region_data) <- seqlengths(tx)
# plan(sequential)

###----------5.3.3 H3K27me3 ----------------------
tmp.files <- list.files(path = "res/chipseq-process/H3K27me3/bedgraph/",full.names = T)
plan(multiprocess, workers = length(tmp.files))
future_lapply(tmp.files,FUN = function(x){
  tmp.samples <- gsub(pattern = "res/chipseq-process/H3K27me3/bedgraph/GSE76687_",replacement = "",x)
  tmp.samples <- gsub(pattern = ".bedgraph",replacement = ".bw",tmp.samples)
  tmp.output <- paste0("res/chipseq-process/data_to_plot/H3K27me3/GSE76687_",tmp.samples)
  myBegraphToBigwig(tmp.input = x,
                    tmp.output = tmp.output,
                    region = NULL,
                    tmp.levels = seqlevels(tx),
                    tmp.seqlength = seqlengths(tx))
})
plan(sequential)



####------5.4 plot track-----------


#####-----5.4.1 make it as a function-------
tmp.track.plot <- function(tmp.region=tmp.region,tmp.ylim=c(0,3),tmp.type="bar",downsample.rate=0.1){
  p_gene <- myGenePlot(annotation = tx,
                       region = tmp.region,
                       arrow_sbreaks = 600,
                       font_size = 18,label_size = 8)+
    xlab(NULL)+
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          axis.title.y = element_text(angle = 0,vjust = 0.5,hjust = 0.5),
          axis.line = element_blank(),axis.text.x = element_blank(),
          axis.ticks = element_blank())
  #p_gene
  
  tmp.dir <- "res/atac-seq-process/data/bigwig/"
  tmp.files <- list.files(path = tmp.dir,pattern = "bw",full.names = T)
  tmp.levels <- c("2cell_early","2cell","8cell","ICM")
  tmp.files <- paste0("res/atac-seq-process/data/bigwig/GSE66581_atac_",tmp.levels,".bw")
  seq_along(tmp.files)
  ###tmp.colors <- shaman.color(length(tmp.files))
  tmp.colors <- c("#4d311a","#d29156","#6d82a3","#6e648a")
  #scales::show_col(tmp.colors)
  
  p_ATAC_track_list <- lapply(seq_along(tmp.files),function(ii){
    cat(ii,sep = "\n")
    p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                       bigwig = tmp.files[ii],
                       smooth = 100,
                       lognorm = T,
                       type = tmp.type,
                       y_label = tmp.levels[ii],
                       fontsize=18,
                       track.color=tmp.colors[ii],
                       tmp.ylimits=tmp.ylim,
                       max.downsample = 3000,
                       downsample.rate = downsample.rate,
                       tmp.seed=42)
    
    return(p)
  })
  
  
  tmp.dir <- "res/chipseq-process/data_to_plot/H3K4me3/"
  tmp.files <- list.files(path = tmp.dir,pattern = "bw",full.names = T)
  tmp.levels <- c("MII_oocyte","PN5_zygote","2cell_early","2cell_late","4cell","8cell","ICM")
  tmp.files <- paste0("res/chipseq-process/data_to_plot/H3K4me3/GSE71434_",tmp.levels,"_K4me3.bw")
  seq_along(tmp.files)
  ###tmp.colors <- shaman.color(length(tmp.files))
  tmp.colors <- c("#b4d1db","#6c1c1f","#4d311a","#d29156","#7f8a65","#6d82a3","#6e648a")
  #scales::show_col(tmp.colors)
  
  
  p_H3K4me3_track_list <- lapply(seq_along(tmp.files),function(ii){
    cat(ii,sep = "\n")
    p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                       bigwig = tmp.files[ii],
                       smooth = 100,
                       lognorm = T,
                       type = tmp.type,
                       y_label = tmp.levels[ii],
                       fontsize=18,
                       track.color=tmp.colors[ii],
                       tmp.ylimits=tmp.ylim,
                       max.downsample = 3000,
                       downsample.rate = downsample.rate,
                       tmp.seed=42)
    return(p)
  })
  
  
  tmp.dir <- "res/chipseq-process/data_to_plot/H3K27me3/"
  tmp.files <- list.files(path = tmp.dir,pattern = "bw",full.names = T)
  tmp.levels <- c("MII_oocyte","PN5_zygote","2cell","ICM","E6.5_epiblast")
  tmp.files <- paste0("res/chipseq-process/data_to_plot/H3K27me3/GSE76687_",tmp.levels,"_k27me3.bw")
  seq_along(tmp.files)
  ###tmp.colors <- shaman.color(length(tmp.files))
  tmp.colors <- c("#b4d1db","#6c1c1f","#d29156","#6e648a","#7aada3")
  #scales::show_col(tmp.colors)
  
  
  p_H3K27me3_track_list <- lapply(seq_along(tmp.files),function(ii){
    cat(ii,sep = "\n")
    p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
                       bigwig = tmp.files[ii],
                       smooth = 100,
                       lognorm = T,
                       type = tmp.type,
                       y_label = tmp.levels[ii],
                       fontsize=18,
                       track.color=tmp.colors[ii],
                       tmp.ylimits=tmp.ylim,
                       max.downsample = 3000,
                       downsample.rate = downsample.rate,
                       tmp.seed=42)
    return(p)
  })
  
  
  p_track_list <- c(p_ATAC_track_list,
                    p_H3K4me3_track_list,
                    p_H3K27me3_track_list)
  
  p <- patchwork::wrap_plots(p_track_list,ncol = 1) +
    plot_annotation(title = tmp.region) & 
    theme(plot.title = element_text(hjust = 0.5,size = 18),
          plot.margin = unit(c(0,0,0,0), "cm"),
          axis.title.y = element_text(size = 18,angle = 0,vjust = 0.5,hjust = -1),
          axis.ticks.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_blank()) 
  p1 <- p / p_gene
  p1
  return(p1)
}

####------5.4.2 plot track--------

tmp.region <- "chr17:35640977-35648722"
tmp.region <- "chr17:35638000-35648722"
p1 <- tmp.track.plot(tmp.region = tmp.region,
                     tmp.ylim = c(0,1),
                     tmp.type = "coverage",
                     downsample.rate = 0.3)

# tss <- c(35640977-1000,35640977+1000)
# p1 & geom_rect(aes(xmin=tss[1],xmax=tss[2],ymin=-Inf,ymax=Inf),fill="yellow")
p1
ggplot() + 
  geom_vline(xintercept = 35640977)+
  scale_x_continuous(limits = c(35640977,35648722))+
  geom_rect(aes(xmin=tss[1],xmax=tss[2],ymin=0,ymax=1,fill="blue"))+
  theme_cowplot()

tss <- c(35640977-2000,35640977-2000)

ggsave(filename = myFileName(prefix = "res/fig/fig1_show_pou5f1_track",suffix = ".png"),
       width = 16,height = 8)
ggsave(filename = myFileName(prefix = "res/fig/fig1_show_pou5f1_track",suffix = ".pdf"),
       width = 16,height = 8)



###-------5.4.3 show RNA-seq----------

####-------5.4.3.1 load data--------
##### RNA-seq xiewei data
xiewei.data.exp <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myNormalize(myRemoveNA(xiewei.data.exp))

beforeEPI.exp.mean <- readRDS(file = "res/R/beforeEPI_exp_mean_20210114.rds")
gene_symbols <- intersect(rownames(xiewei.data.exp),rownames(beforeEPI.exp.mean))

####rename for plot
colnames(beforeEPI.exp.mean)[1] <- "MII_oocyte"
colnames(xiewei.data.exp)[4] <- "late2cell"
colnames(xiewei.data.exp)[3] <- "early2cell"
print(colnames(xiewei.data.exp),quote=T)
beforeEPI.exp.mean$ICM <- beforeEPI.exp.mean %>%
  dplyr::select(c("earlyblast","midblast","lateblast")) %>%
  rowMeans()


####load Epiblast data

seu <- readRDS("res/R/early_gastrulation_20210115.rds")
seu@assays$RNA@counts <- 2^(seu@assays$RNA@data)-1

levels(seu) <- rev(c("EPI","ExE","VE"))
TSNEPlot(seu,label=T)

colnames(seu[[]])
tmp.seu <- subset(seu,cell_type=="EPI" & Stage == "E6.5")
E6.5_EPI.mean <- AverageExpression(tmp.seu,slot = "counts")$RNA

tmp.exp <- log2(E6.5_EPI.mean+1)
tmp.exp <- tmp.exp[gene_symbols,]


####----5.4.3.2 prepare data---------
Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
data.plot.1 <- beforeEPI.exp.mean[gene_symbols,Stage]
data.plot.2 <- xiewei.data.exp[gene_symbols,Stage]
data.plot.1$E6.5.EPI <- tmp.exp
tmp.gene <- "Pou5f1"

boxplot(data.plot.1)
boxplot(data.plot.2)

tmp.data.plot.1 <- data.frame(value=t(data.plot.1[tmp.gene,]),
                              type="scRNA-seq",
                              gene=tmp.gene,stringsAsFactors = F) %>%
  rownames_to_column("stage")

tmp.data.plot.2 <- data.frame(value=t(data.plot.2[tmp.gene,]),
                              type="RNA-seq",
                              gene=tmp.gene,stringsAsFactors = F) %>%
  rownames_to_column("stage")
data.plot <- rbind(tmp.data.plot.1,tmp.data.plot.2)
colnames(data.plot)[2] <- "value"
data.plot <- rbind(data.plot,
                   data.frame(stage="E6.5.EPI",
                              value=NA,
                              type="RNA-seq",
                              gene="Pou5f1"))

#####-----5.4.4.4 let's plot data----------


plot.list <- list()
#### Match ATAC
tmp.stage <- c("early2cell","late2cell","8cell","ICM")
tmp.colors <- c("#4d311a","#d29156","#6d82a3","#6e648a")
names(tmp.colors) <- tmp.stage
tmp.plot.list <- lapply(tmp.stage,FUN = function(ii){
  tmp.data.plot <- data.plot %>%
    dplyr::filter(stage == ii)
  tmp.data.plot
  p <- ggplot(tmp.data.plot,aes(type,value))+
    geom_bar(stat = "identity",fill=tmp.colors[ii])+
    xlab(NULL)+
    ylab("log2(RPKM+1)")+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,10),
                       position = "right")+
    theme_cowplot(font_size = 18)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(p)
})
plot.list <- c(plot.list,tmp.plot.list)

####Match K4Me3

tmp.stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
tmp.colors <- tmp.colors <- c("#b4d1db","#6c1c1f","#4d311a","#d29156","#7f8a65","#6d82a3","#6e648a")
names(tmp.colors) <- tmp.stage
tmp.plot.list <- lapply(tmp.stage,FUN = function(ii){
  tmp.data.plot <- data.plot %>%
    dplyr::filter(stage == ii)
  tmp.data.plot
  p <- ggplot(tmp.data.plot,aes(type,value))+
    geom_bar(stat = "identity",fill=tmp.colors[ii])+
    xlab(NULL)+
    ylab("log2(RPKM+1)")+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,10),
                       position = "right")+
    theme_cowplot(font_size = 18)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(p)
})
plot.list <- c(plot.list,tmp.plot.list)

####Match K27me3

tmp.stage <- c("MII_oocyte","zygote","late2cell","ICM","E6.5.EPI")
tmp.colors <-  c("#b4d1db","#6c1c1f","#d29156","#6e648a","#7aada3")
names(tmp.colors) <- tmp.stage
tmp.plot.list <- lapply(tmp.stage,FUN = function(ii){
  tmp.data.plot <- data.plot %>%
    dplyr::filter(stage == ii)
  tmp.data.plot
  p <- ggplot(tmp.data.plot,aes(type,value))+
    geom_bar(stat = "identity",fill=tmp.colors[ii])+
    xlab(NULL)+
    ylab("log2(RPKM+1)")+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,10),
                       position = "right")+
    theme_cowplot(font_size = 18)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(p)
})
plot.list <- c(plot.list,tmp.plot.list)
p <- wrap_plots(plot.list,ncol = 1) & 
  ylab(NULL) &
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,10),
                     breaks = seq(0,10,length.out=3),
                     position = "right") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = myFileName(prefix = "res/fig/fig1_show_pou5f1_RNA-seq_track",
                             suffix = ".png"),
       width = 1.5,height = 8)
ggsave(filename = myFileName(prefix = "res/fig/fig1_show_pou5f1_RNA-seq_track",
                             suffix = ".pdf"),
       width = 1.5,height = 8)



