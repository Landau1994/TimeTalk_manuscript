####time-talk
###time-talk 3D fig4;
####how autocrine influence cell signaling?
####20210420
####how genome structure influence cell signaling?
####how signaling influence genome structure?
#### using Promoter as gene coordinate
#### make it clear
####author: wlt
####20210406
####20210419
####20210505
####20210510, Give up; or keep
####20210517， work hard for validate
####20210531, step by step configuration;
####20210712, keep going, go foward;
####20210719, should I give up or just keep chasing pavements

####20210809, I hope finish this part in this week
####20210812, from ATAC-seq, K4, et.c to illustrate Add epigenomics data
####20210813, Don't think Hi-C data, regulatory network

####20211101, The past is writen, but the future is left for us to write
####----------0. load requirement environment----------------

####---------0.1 load package-------------------------------
library(tidyverse)
library(Seurat)
library(ggplot2)
library(HiTC)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(cowplot)
library(VennDiagram)
library(ggbeeswarm)
library(ggsignif)
library(ggbio)
library(HiTC)
library(GenomicFeatures)
library(RcisTarget)
library(ggrepel)
library(fields)
library(ggthemes)
library(msigdbr)
library(org.Mm.eg.db)
library(CellChat)
library(future.apply)
library(UpSetR)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(ComplexHeatmap)
source("code/myUtils.R")
options(connectionObserver = NULL)


####----0.2 define color------
rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))) 
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
col.spectral <- colorRampPalette(brewer.pal(11,'Spectral')[-6])
test.color.3 <- colorRampPalette(c("#f86e11","#e9970a","#71a701","#62b474","#71c3ac","#9fc4ca"))
rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))
skyblueYelow <- colorRampPalette(c("skyblue","black","yellow"))
skybluered <- colorRampPalette(c("skyblue","black","orange"))
solarExtra <- colorRampPalette(c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
hic.red <- colorRampPalette(c("white","red"))
hic.orage <- colorRampPalette(brewer.pal(9,"Oranges"))
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
coldwarm <- colorRampPalette(colors = mypalette)

bkcolor <- c(colorRampPalette(c(brewer.pal(9,"Blues" )[4:9],"#1a1919"))(50),
             colorRampPalette(c("#1a1919",rev(brewer.pal( 9,"YlOrBr" ))[1:6]))(50))

bkcolor <- colorRampPalette(colors = bkcolor)
mycolor.bar(bkcolor(100),min = -1)

#####--------1. plot Is heat map-------------


####---------1.1 load eLR-----------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
cutoff <- 0.1
tmp.res.LR.pos <- tmp.res.cor.1 %>%
  dplyr::filter(PCC > cutoff)
tmp.res.LR.neg <- tmp.res.cor.1 %>%
  dplyr::filter(PCC < -cutoff)
tmp.res.LR <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff)


####load mat
IS.mat <- readRDS(file = "res/R/early.scRNAseq_cytotalk_LR_IS_20210205.rds")

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5")
####LR.pairs
tmp.mat.plot <- IS.mat[tmp.res.LR$LRpairs,Stage.levels]
tmp.mat.plot <- log10(tmp.mat.plot+1)
tmp.mat.plot <- myRemoveNA(tmp.mat.plot)


head(tmp.mat.plot)

#####------1.2 mfuzz analysis---------
library(Mfuzz)
tmp.mat.plot <- as.matrix(tmp.mat.plot)
eset <- new("ExpressionSet",exprs = tmp.mat.plot)

### filter gene 
gene.f <- fill.NA(eset,mode="mean")
tmp <- filter.std(gene.f,min.std=0)
### standard
gene.s <- standardise(tmp)

set.seed(4)
c <- 6
m <- mestimate(gene.s)
## clustering
cl <- mfuzz(gene.s,centers = c, m = m)
# ## look cluster
# mfuzz.plot(eset = gene.s,cl)

## show ID
cl$cluster[cl$cluster == 1]
##  OUTPUT ID
###write.table(cl$cluster,"output.txt",quote=F,row.names=T,col.names=F,sep="\t")
## plot curve



ttttt <- gene.s@assayData$exprs
cl$cluster[cl$cluster == 1]

annotation_row <- as.data.frame(cl$cluster)
colnames(annotation_row) <- "cluster"
annotation_row$cluster <- paste0("cluster_",annotation_row$cluster)
annotation_row <- annotation_row %>%
  arrange(cluster)
annotation_row$cluster <- plyr::mapvalues(annotation_row$cluster,
                                          from = paste0("cluster_",1:6),
                                          to = paste0("cluster_",c(1,6,5,4,2,3))) 
annotation_row <- annotation_row %>%
  arrange(cluster)

saveRDS(object = annotation_row,
        file = myFileName(prefix = "res/R/eLR_IS_timecourse_mfuzz_cluster",suffix = ".rds"))
saveRDS(object = ttttt,
        file = myFileName(prefix = "res/R/eLR_IS_mfuzz_scale_mat",suffix = ".rds"))


#####------1.3 load data-------
tmp.cluster.color <- ggsci::pal_npg()(6)
scales::show_col(tmp.cluster.color)
names(tmp.cluster.color) <- paste0("cluster_",1:6)
ann_colors <- list(cluster = tmp.cluster.color)  
annotation_row <- readRDS(file = "res/R/eLR_IS_timecourse_mfuzz_cluster_20210611.rds")
ttttt <- readRDS(file = "res/R/eLR_IS_mfuzz_scale_mat_20210611.rds")



ph2 <- ComplexHeatmap::pheatmap(ttttt[rownames(annotation_row),],
                                cluster_cols = F,
                                cluster_rows = F,
                                color = rdbu(100),
                                show_rownames = F,
                                annotation_colors = ann_colors,
                                annotation_row = annotation_row,
                                annotation_names_row = F,
                                angle_col = "45",
                                name = "zscore",
                                fontsize = 20,
                                heatmap_legend_param = list(title_gp=gpar(fontsize=20,fontface="bold")),
                                labels_row = NULL,
                                border_color = NA)

png(filename = myFileName(prefix = "res/fig/fig4_IS_at_different_stage_rdbu",suffix = ".png"),
    width = 6,height = 8,units = "in",res = 350,bg = "white")
ph2
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig4_IS_at_different_stage_rdbu",suffix = ".pdf"),
    width = 6,height = 8)
ph2
dev.off()


####-------1.4 plot time curve----------------
tmp.clusters <- unique(annotation_row$cluster)
plot.list <- list()
for(ii in 1:6){
  cat(tmp.clusters[ii],sep = "\n")
  tmp.id <- annotation_row %>%
    dplyr::filter(cluster == tmp.clusters[ii]) %>%
    rownames()
  
  
  data.plot <- ttttt[tmp.id,] %>%
    as.data.frame() %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "Stage",value="zscore",-LRpairs) %>%
    mutate(Stage=factor(Stage,levels = Stage.levels))
  
  
  tmp.color = ggsci::pal_npg()(6)
  range(data.plot$zscore)
  
  p <- ggplot(data.plot,aes(Stage,zscore))+
    geom_point(alpha=0.03,color=tmp.color[ii])+
    geom_line(aes(group=LRpairs),color=tmp.color[ii],alpha=0.3)+
    geom_smooth(method = "gam",aes(group=""),se = T,color="black")+
    xlab(NULL)+
    scale_y_continuous(limits = c(-3,4))+
    ggtitle(tmp.clusters[ii])+
    theme_cowplot(font_size = 22)+
    theme(legend.position = "none",
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust = 1))
  plot.list <- c(plot.list,list(p))
}
plot_grid(plotlist = plot.list,nrow = 2)
ggsave(filename = myFileName(prefix = "res/fig/fig4_eLR_cluster_expression",suffix = ".png"),
       width = 16,height = 9,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/fig4_eLR_cluster_expression",suffix = ".pdf"),
       width = 16,height = 9,dpi = 350)

####-----1.5 PCC -----------

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
cutoff <- 0.1
tmp.res.LR <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff)
head(tmp.res.LR)

tmp.tttt <- annotation_row %>%
  rownames_to_column("LRpairs")
tmp.data.plot <- merge(tmp.res.LR,tmp.tttt,by="LRpairs")
tmp.data.plot$cluster <- factor(tmp.data.plot$cluster,levels = paste0("cluster_",1:6))

ggplot(tmp.data.plot,aes(cluster,PCC,fill=cluster))+
  geom_boxplot(size=1)+
  scale_fill_manual(values = ggsci::pal_npg()(6))+
  theme_cowplot(font_size = 28)+
  xlab(NULL)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,vjust = 0.5),
        axis.line = element_line(size=2),
        axis.ticks = element_line(size = 2))
ggsave(filename = myFileName(prefix = "res/fig/fig4_eLR_PCC",suffix = ".png"),
       width = 8,height = 8,dpi = 350)

####--------2.maternal factor------
####----2.1.1 load eLR data -----
annotation_row <- readRDS(file = "res/R/eLR_IS_timecourse_mfuzz_cluster_20210611.rds")
eLR.df <- annotation_row %>%
  rownames_to_column("eLR")
####----2.1.2 load maternal factor------
tmp.df <- readxl::read_excel(path = "database/DBTMEE/maternal_gene.xlsx",skip = 2) %>%
  dplyr::filter(Cluster == "Maternal RNA")
maternal.gene <- tmp.df$Gene
####---2.1.3 load ZGA gene-------
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)



####-----2.2 maternal overlap---------
####----2.2.1 overlap with maternal gene---------
tmp.list <- lapply(unique(eLR.df$cluster),FUN = function(ii){
  tmp.x <- eLR.df %>%
    dplyr::filter(cluster==ii) %>%
    pull(eLR)
  Lgenelist <- unlist(lapply(strsplit(tmp.x,split = "_"),FUN = function(x) x[1]))
  Rgenelist <- unlist(lapply(strsplit(tmp.x,split = "_"),FUN = function(x) x[2]))
  Lgene <- unique(Lgenelist)
  Rgene <- unique(Rgenelist)
  LRgene <- c(Lgene,Rgene)
  tmp.res <- intersect(maternal.gene,LRgene)
  tmp.res.ratio <- length(tmp.res) / length(LRgene)
  return(tmp.res.ratio)
})
names(tmp.list) <- unique(eLR.df$cluster)
tmp.list <- unlist(tmp.list)
### plot 
png(filename = myFileName(prefix = "res/fig/fig5_maternal_ratio",suffix = ".png"),
    width = 8,height = 8,units = "in",res = 350)
mp <- barplot(tmp.list,col = tmp.cluster.color,
              cex.axis=1.5,
              xaxt="n",
              srt=45,
              cex.lab=1.5,
              font.lab = 2,
              ylab = "ratio",
              main = "Maternal factor overlap ratio",
              cex.main=2,
              lwd=2)
tot <- tmp.list
text(mp, tot+0.005, format(tot,digits = 2), 
     xpd = TRUE, 
     col = "black",
     cex=1.5)
text(mp, par("usr")[3], 
     labels = names(tmp.list), 
     srt = 45, 
     adj = c(1.0,1.0), 
     xpd = TRUE, 
     cex=1.5)
dev.off()

#####------2.2.2 overlap with ZGA gene---------

tmp.list <- lapply(unique(eLR.df$cluster),FUN = function(ii){
  tmp.x <- eLR.df %>%
    dplyr::filter(cluster==ii) %>%
    pull(eLR)
  Lgenelist <- unlist(lapply(strsplit(tmp.x,split = "_"),FUN = function(x) x[1]))
  Rgenelist <- unlist(lapply(strsplit(tmp.x,split = "_"),FUN = function(x) x[2]))
  Lgene <- unique(Lgenelist)
  Rgene <- unique(Rgenelist)
  LRgene <- c(Lgene,Rgene)
  tmp.res <- intersect(zga_gene_zhangyi,LRgene)
  tmp.res.ratio <- length(tmp.res) / length(LRgene)
  return(tmp.res.ratio)
})
names(tmp.list) <- unique(eLR.df$cluster)
tmp.list <- unlist(tmp.list)
### plot 
png(filename = myFileName(prefix = "res/fig/fig5_ZGA_ratio",suffix = ".png"),
    width = 8,height = 8,units = "in",res = 350)
mp <- barplot(tmp.list,col = tmp.cluster.color,
              cex.axis=1.5,
              xaxt="n",
              srt=45,
              cex.lab=1.5,
              font.lab = 2,
              ylab = "ratio",
              main = "ZGA gene overlap ratio",
              cex.main=2,
              lwd=2)
tot <- tmp.list
text(mp, tot+0.005, format(tot,digits = 2), 
     xpd = TRUE, 
     col = "black",
     cex=1.5)
text(mp, par("usr")[3], 
     labels = names(tmp.list), 
     srt = 45, 
     adj = c(1.0,1.0), 
     xpd = TRUE, 
     cex=1.5)
dev.off()


#####----overlap of essential and housekeeping gene----------
tmp.dir <- "database/deg_annotation_e.csv/deg_annotation_e.csv"
tmp.df <- read.delim(file = tmp.dir,sep = ";",stringsAsFactors = F,header = F)
essential.gene <- tmp.df %>%
  dplyr::filter(V8 == "Mus musculus") %>%
  pull(V3)
### load HK gene
tmp.dir <- "database/HRT_housekeeping_gene_database/Housekeeping_GenesMouse.csv"
tmp.df <- read.delim(file = tmp.dir,sep = ";",stringsAsFactors = F,header = T)
HK.gene <- tmp.df %>%
  pull(Gene) %>%
  unique()

tmp.list <- list(essential.gene=essential.gene,
                 HK.gene=HK.gene)
p <- venn.diagram(tmp.list,filename = NULL)
grid.draw(p)





