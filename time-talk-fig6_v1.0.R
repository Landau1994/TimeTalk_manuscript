####Fig6
####20210922
####count Ligand receptor gene ratio
####20211124, revised version

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
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
library(RcisTarget)
library(ggprism)
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
tmp.color <- rev(c("#f4ea0c","#ffc807",
                   "#f6901f","#f05e23",
                   "#ed3324","#ec1c25",
                   "#bc243c","#8d3966",
                   "#5a5391","#366cb3",
                   "#1f65af","#174c81",
                   "#163256","#0e192a",
                   "#000000"))
bkro <- colorRampPalette(colors = tmp.color)
mycolor.bar(bkro(100),min = -1)
shaman.color <- c(c("#08306b","#084d96","#3f51b5","#2196f3","#03a9f4","#00bcd4"),
                  colorRampPalette(c("#00bcd4","lightgrey","#ffeb3b"))(5)[2:4],
                  rev(c("#800026","#B60026","#ff5722","#ff9800","#ffc107","#ffeb3b")))

shaman.color <- colorRampPalette(colors = shaman.color)

molandi.color <- c("#e4d4c5","#c6b1ac","#764e56",
                   "#f9d9c0","#d19477","#93675a",
                   "#f0e9cd","#b9a783","#796656",
                   "#cdc1b1","#a2967e","#656356",
                   "#d8e7e4","#9eb2b1","#5a6873",
                   "#ccd8b0","#7e8563","#50463d")
scales::show_col(molandi.color,ncol = 3)

molandi.color <- colorRampPalette(colors = molandi.color)
scales::show_col(molandi.color(36),ncol = 6)

pheatmap::pheatmap(farver::decode_colour(shaman.color(36), "rgb", "hcl"),
                   cluster_rows = F,cluster_cols = F)


####-----------1.Overexpression KLF4 RNA-seq data-------------


tmp.files <- "data/GSE90894_RPKM_mRNAseq_table.txt"
tmp.df <- read.delim(file = tmp.files,stringsAsFactors = F,skip = 4)
tmp.idx <- c("name","MEFs","MEFs_Klf4")
tmp.df <- tmp.df[,tmp.idx]
tmp.fc.df <- tmp.df %>%
  mutate(log2fc=log2(MEFs_Klf4/MEFs))
mm9.gene <- tmp.df$name

tmp.df <- read.delim(file = "res/txt/eLRgene_network_20210512.txt",stringsAsFactors = F,header = T)
head(tmp.df)
tmp.gene <- tmp.df %>%
  dplyr::filter(source == "Klf4") %>%
  pull(target) %>%
  unique()

tmp.fc.df.plot <- tmp.fc.df %>%
  dplyr::filter(name %in% tmp.gene) %>%
  arrange(-log2fc)
tmp.max <- max(range(tmp.fc.df.plot$log2fc,finite=T))

tmp.fc.df.plot$log2fc[is.na(tmp.fc.df.plot$log2fc)] <- 0
tmp.fc.df.plot$log2fc[tmp.fc.df.plot$log2fc==Inf] <- tmp.max
tmp.fc.df.plot$log2fc[tmp.fc.df.plot$log2fc==-Inf] <- -tmp.max



tmp.fc.df.plot <- tmp.fc.df.plot %>%
  arrange(-log2fc)
tmp.fc.df.plot$group <- ifelse(tmp.fc.df.plot$log2fc>0,"up","down")
tmp.fc.df.plot$group <- factor(tmp.fc.df.plot$group,levels = c("up","down"))
tmp.fc.df.plot$name <- factor(tmp.fc.df.plot$name,
                              levels = tmp.fc.df.plot$name)


ggplot(tmp.fc.df.plot,aes(name,log2fc,fill=group))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#b2182b","#2166ac"))+
  ylab("log2(fold change)")+
  xlab("gene symbol")+
  theme_few(base_size = 16)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90,vjust = 0.5))
ggsave(filename = myFileName(prefix = "res/fig/KLF4_OE_foldchange",suffix = ".jpg"),
       width = 18,height = 8,dpi = 350)



data.df.2 <- tmp.fc.df.plot[,c("name","log2fc")]

####----------2.KO data----------

###ESC KLF4 KO, GSE129495

###load package
library(DESeq2)
library(EnsDb.Mmusculus.v79)
###prepare gene track
mm10_gene <- genes(EnsDb.Mmusculus.v79)
mm10_gene <- as.data.frame(mm10_gene)
mm10_gene <- mm10_gene %>%
  dplyr::filter(gene_biotype == "protein_coding")

###load counts
data.df <- read.delim(file = "data/GSE129495_raw.tsv.gz",stringsAsFactors = F)
data.df <- data.df %>% 
  column_to_rownames("Ens_ID")
tmp.gene <- intersect(rownames(data.df),mm10_gene$gene_id)
data.df <- data.df[tmp.gene,]

###meta data
data.meta <- data.frame(sample_id = colnames(data.df),condition=c("WT","WT","KO","KO"))
data.meta <- data.meta %>%
  column_to_rownames("sample_id")
data.meta$condition <- factor(data.meta$condition,levels = c("WT","KO"))

####bulid cds
cds <- DESeqDataSetFromMatrix(countData = data.df,
                              colData = data.meta,
                              design = ~ condition)

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

####add ranges
Ens.ghs <- genes(EnsDb.Mmusculus.v79)
Ens.ghs.subset <- Ens.ghs[tmp.gene,]
rowRanges(cds) <- Ens.ghs.subset  
cds <- DESeq(cds)

####get fpkm
exprs.fpkm <- fpkm(cds)
rownames(exprs.fpkm) <- plyr::mapvalues(rownames(exprs.fpkm),
                              from = mm10_gene$gene_id,
                              to = mm10_gene$symbol,
                              warn_missing = F)

####get degene
head(exprs.fpkm[,1:3])
tmp.res.diff <- as.data.frame(results(cds)) %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(-log2FoldChange)
data.df <- tmp.res.diff
rownames(data.df) <- plyr::mapvalues(rownames(data.df),
                                     from = mm10_gene$gene_id,
                                     to = mm10_gene$symbol,
                                     warn_missing = F)

# ####------2.1 overlap with DEgene----------
# 
# tmp.df <- read.delim(file = "res/txt/eLRgene_network_20210512.txt",stringsAsFactors = F,header = T)
# head(tmp.df)
# tmp.gene <- tmp.df %>%
#   dplyr::filter(source == "Klf4") %>%
#   pull(target) %>%
#   unique()
# 
# #tmp.gene <- sample(mm9.gene,length(tmp.gene))
# tmp.overlap <- intersect(tmp.gene,data.df$symbol)
# tmp.ratio <- length(tmp.overlap)/length(data.df$symbol)
# 
# 
# 
# plan(multiprocess, workers = 10)
# tmp.list <- future_replicate(n = 1000,expr = {
#   tmp.gene <- sample(mm9.gene,length(tmp.gene))
#   tmp.overlap <- intersect(tmp.gene,data.df$symbol)
#   tmp.ratio <- length(tmp.overlap)/length(data.df$symbol)
#   return(tmp.ratio)
# })
# 
# tmp.p <- length(tmp.list[tmp.list > tmp.ratio]) / length(tmp.list)
# 
# 
# 
# tmp.list <- list(Klf4_target=tmp.gene,DE_gene=data.df$symbol)
# 
# png(filename = myFileName(prefix = "res/fig/Fig6_Klf4_KO_",suffix = ".png"),width = 8,height = 8,units = "in",bg = "white",res = 350)
# p <- VennDiagram::venn.diagram(tmp.list,
#                                fill = c("blue", "darkred"),
#                                filename = NULL,
#                                total.population = length(mm9.gene),
#                                hyper.test = T,
#                                main="Klf4 KO",main.fontface = "bold",
#                                lower.tail = F)
# 
# grid::grid.draw(p)
# dev.off()

####--------2.2 LR ratio-------------

####load LR gene 
#####------2.2.1 LRgene -------
mm10_gene.set <- mm10_gene$symbol

LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm10_gene.set)

DE_gene <- rownames(data.df)

LRgene.1 <- LRgene

####phyper(k-1,M,N-M,n)
k <- length(intersect(LRgene,DE_gene))
M <- length(DE_gene)
N <- length(mm10_gene.set)
n <- length(LRgene)
x <- M/N
y <- k/n

pvalue <- phyper(k-1,M,N-M,n,lower.tail = F)
tmp.df <- data.frame(cat = c("Total","LRgene"),
                 ratio = c(x,y),
                 stringsAsFactors = F)
data.plot <- tmp.df
pvalue.list <- NULL
pvalue.list <- c(pvalue.list,pvalue)

######---------2.2.2 eLR ratio-----------
tmp.df <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
LRpairs <- tmp.df %>%
  dplyr::filter(abs(PCC) > 0.1) %>%
  pull(LRpairs)

Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm10_gene.set)

eLR <- LRgene
DE_gene <- rownames(data.df)
mm10_gene.set <- mm10_gene$symbol


####phyper(k-1,M,N-M,n)
k <- length(intersect(LRgene,DE_gene))
M <- length(DE_gene)
N <- length(mm10_gene.set)
n <- length(LRgene)

x <- M/N
y <- k/n

pvalue <- phyper(k-1,M,N-M,n,lower.tail = F)


tmp.df <- data.frame(cat = c("eLR"),
                        ratio = y,
                        stringsAsFactors = F)
data.plot <- rbind(data.plot,tmp.df)
pvalue.list <- c(pvalue.list,pvalue)

####-----------2.2.3 Klf4 target eLR ratio--------------
tmp.df <- read.delim(file = "res/txt/eLRgene_network_20210512.txt",stringsAsFactors = F,header = T)
head(tmp.df)
tmp.gene <- tmp.df %>%
  dplyr::filter(source == "Klf4") %>%
  pull(target) %>%
  unique()
mm10_gene.set <- mm10_gene$symbol
LRgene <- intersect(tmp.gene,mm10_gene.set)

Klf4_eLR_target <- LRgene
DE_gene <- rownames(data.df)

####phyper(k-1,M,N-M,n)
k <- length(intersect(LRgene,DE_gene))
M <- length(DE_gene)
N <- length(mm10_gene.set)
n <- length(LRgene)

x <- M/N
y <- k/n


pvalue <- phyper(k-1,M,N-M,n,lower.tail = F)


tmp.df <- data.frame(cat = c("Klf4_target_eLR"),
                        ratio = y,
                        stringsAsFactors = F)
data.plot <- rbind(data.plot,tmp.df)
pvalue.list <- c(pvalue.list,pvalue)


data.plot$cat <- factor(data.plot$cat,levels = c("Total","LRgene","eLR","Klf4_target_eLR"))
####------2.2.4 bar plot---------

ggplot(data = data.plot,aes(cat,ratio,fill=cat))+
  geom_text(aes(label=round(ratio,2),vjust=-0.5),size=12)+
  xlab(NULL)+
  ggtitle("Klf4 KO DEgene \n in LRgene enrichment")+
  geom_bar(stat="identity")+
  geom_signif(comparisons = list(c("Total","LRgene"),
                                 c("Total","eLR"),
                                 c("Total","Klf4_target_eLR")),
              annotations = paste0("p=",signif(pvalue.list,3)),
              y_position = c(0.15,0.22,0.25)+0.05,
              size = 2,
              tip_length = 0.1,
              textsize = 10)+
  scale_fill_manual(values = c("grey",blues(4)[-1]))+
  scale_x_discrete(labels = c("all","LR","eLR","Klf4 target \n eLR"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.4),breaks = seq(0,0.4,0.1))+
  theme_cowplot(font_size = 28)+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

ggsave(filename = myFileName(prefix = "res/fig/Fig6_Klf4_KO_DEgene_in_LRgene_value",suffix = ".png"),
       width = 8,height = 8,dpi = 350,units = "in",bg = "white")

### p2star version
# ggplot(data = data.plot,aes(cat,ratio,fill=cat))+
#   geom_text(aes(label=round(ratio,2),vjust=-0.5),size=6)+
#   xlab(NULL)+
#   ggtitle("Klf4 KO DEgene in LRgene enrichment")+
#   geom_bar(stat="identity")+
#   geom_signif(comparisons = list(c("Total","LRgene"),
#                                  c("Total","eLR"),
#                                  c("Total","Klf4_target_eLR")),
#               annotations = as.character(p2star(pvalue.list)),
#               y_position = c(0.15,0.22,0.25),
#               size = 1,
#               textsize = 6)+
#   scale_fill_manual(values = c("grey",blues(4)[-1]))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,0.3),breaks = seq(0,0.3,0.1))+
#   theme_cowplot(font_size = 18)+
#   theme(axis.line = element_line(size = 1),
#         plot.title = element_text(hjust = 0.5),
#         legend.position="none")
# 
# ggsave(filename = myFileName(prefix = "res/fig/Fig6_Klf4_KO_DEgene_in_LRgene",suffix = ".png"),
#        width = 8,height = 8,dpi = 350,units = "in",bg = "white")


#####-------2.2.5 venndigram-------

tmp.list <- list(DE = DE_gene,
                 LR = LRgene.1,
                 eLR = eLR,
                 Klf4_target_eLR = Klf4_eLR_target)
png(filename = myFileName(prefix = "res/fig/Fig6_Klf4_KO_DE_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)

plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = T,
     label=T,
     fill=c("white",blues(4)[-1]),
     col="black",
     bg="grey")
dev.off()


####----2.2.6 barplot revised------


pvalue.list <- NULL
pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(DE_gene,mm10_gene.set)),
                                                  N = length(mm10_gene.set),
                                                  k = length(intersect(LRgene.1,DE_gene)),
                                                  n = length(LRgene.1),
                                                  lower.tail = F))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.1,DE_gene)),
                                                  N = length(LRgene.1),
                                                  k = length(intersect(eLR,DE_gene)),
                                                  n = length(eLR),
                                                  lower.tail = F))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(eLR,DE_gene)),
                                                  N = length(eLR),
                                                  k = length(intersect(Klf4_eLR_target,DE_gene)),
                                                  n = length(Klf4_eLR_target),
                                                  lower.tail = F))


pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(DE_gene,mm10_gene.set)),
                                                  N = length(mm10_gene.set),
                                                  k = length(intersect(eLR,DE_gene)),
                                                  n = length(eLR),
                                                  lower.tail = F))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(DE_gene,mm10_gene.set)),
                                                  N = length(mm10_gene.set),
                                                  k = length(intersect(Klf4_eLR_target,DE_gene)),
                                                  n = length(Klf4_eLR_target),
                                                  lower.tail = F))


ggplot(data = data.plot,aes(cat,ratio,fill=cat))+
  geom_text(aes(label=round(ratio,2),vjust=-0.5),size=12)+
  xlab(NULL)+
  ggtitle("KLFs TKO DEgene \n in LRgene enrichment")+
  geom_bar(stat="identity")+
  geom_signif(comparisons = list(c("Total","LRgene"),
                                 c("LRgene","eLR"),
                                 c("eLR","Klf4_target_eLR"),
                                 c("Total","eLR"),
                                 c("Total","Klf4_target_eLR")),
              annotations = paste0("p=",signif(pvalue.list,3)),
              y_position = c(0.20,0.25,0.30,0.35,0.40),
              size = 2,
              textsize = 10)+
  scale_fill_manual(values = c("grey",blues(4)[-1]))+
  scale_x_discrete(labels = c("all","LR","eLR","KLF4 \n target eLR"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.45),breaks = seq(0,0.45,0.1))+
  theme_cowplot(font_size = 28)+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")
ggsave(filename = myFileName(prefix = "res/fig/Fig6_Klf4_DE_DEgene_in_LRgene_value",suffix = ".png"),
       width = 8,height = 8,dpi = 350,units = "in",bg = "white")


####


####------3.OE data---------

####------3.1.load OE data-------------
tmp.files <- "data/GSE90894_RPKM_mRNAseq_table.txt"
tmp.df <- read.delim(file = tmp.files,stringsAsFactors = F,skip = 4)
tmp.idx <- c("name","MEFs","MEFs_Klf4")
tmp.df <- tmp.df[,tmp.idx]
tmp.fc.df <- tmp.df %>%
  mutate(log2fc=log2(MEFs_Klf4/MEFs))
DE_gene <- tmp.fc.df %>%
  dplyr::filter(abs(log2fc) > 1.5) %>%
  pull(name)

DE_gene <- intersect(DE_gene,mm10_gene.set)

####-----3.2 compute ratio-----------
####load LR gene 
#####------3.2.1 LRgene -------
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm10_gene.set)


LRgene.1 <- LRgene

####phyper(k-1,M,N-M,n)
k <- length(intersect(LRgene,DE_gene))
M <- length(DE_gene)
N <- length(mm10_gene.set)
n <- length(LRgene)
x <- M/N
y <- k/n

pvalue <- phyper(k-1,M,N-M,n,lower.tail = F)
tmp.df <- data.frame(cat = c("Total","LRgene"),
                     ratio = c(x,y),
                     stringsAsFactors = F)
data.plot <- tmp.df
pvalue.list <- NULL
pvalue.list <- c(pvalue.list,pvalue)

######---------3.2.2 eLR ratio-----------
tmp.df <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
LRpairs <- tmp.df %>%
  dplyr::filter(abs(PCC) > 0.1) %>%
  pull(LRpairs)

Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm10_gene.set)

eLR <- LRgene


####phyper(k-1,M,N-M,n)
k <- length(intersect(LRgene,DE_gene))
M <- length(DE_gene)
N <- length(mm10_gene.set)
n <- length(LRgene)

x <- M/N
y <- k/n

pvalue <- phyper(k-1,M,N-M,n,lower.tail = F)


tmp.df <- data.frame(cat = c("eLR"),
                     ratio = y,
                     stringsAsFactors = F)
data.plot <- rbind(data.plot,tmp.df)
pvalue.list <- c(pvalue.list,pvalue)

####-----------3.2.3 Klf4 target eLR ratio--------------
tmp.df <- read.delim(file = "res/txt/eLRgene_network_20210512.txt",stringsAsFactors = F,header = T)
head(tmp.df)
tmp.gene <- tmp.df %>%
  dplyr::filter(source == "Klf4") %>%
  pull(target) %>%
  unique()
LRgene <- intersect(tmp.gene,mm10_gene.set)

Klf4_eLR_target <- LRgene

####phyper(k-1,M,N-M,n)
k <- length(intersect(LRgene,DE_gene))
M <- length(DE_gene)
N <- length(mm10_gene.set)
n <- length(LRgene)

x <- M/N
y <- k/n


pvalue <- phyper(k-1,M,N-M,n,lower.tail = F)


tmp.df <- data.frame(cat = c("Klf4_target_eLR"),
                     ratio = y,
                     stringsAsFactors = F)
data.plot <- rbind(data.plot,tmp.df)
pvalue.list <- c(pvalue.list,pvalue)


data.plot$cat <- factor(data.plot$cat,levels = c("Total","LRgene","eLR","Klf4_target_eLR"))

####-----3.3 barplot----
ggplot(data = data.plot,aes(cat,ratio,fill=cat))+
  geom_text(aes(label=round(ratio,2),vjust=-0.5),size=6)+
  xlab(NULL)+
  ggtitle("Klf4 OE DEgene in LRgene enrichment")+
  geom_bar(stat="identity")+
  geom_signif(comparisons = list(c("Total","LRgene"),
                                 c("Total","eLR"),
                                 c("Total","Klf4_target_eLR")),
              annotations = paste0("p=",signif(pvalue.list,3)),
              y_position = c(0.55,0.65,0.75),
              size = 1,
              textsize = 6)+
  scale_fill_manual(values = c("grey",blues(4)[-1]))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.8),breaks = seq(0,0.8,0.1))+
  theme_cowplot(font_size = 28)+
  theme(axis.line = element_line(size = 2),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")
ggsave(filename = myFileName(prefix = "res/fig/Fig6_Klf4_KO_DEgene_in_LRgene_value",suffix = ".png"),
       width = 8,height = 8,dpi = 350,units = "in",bg = "white")

####p2star
ggplot(data = data.plot,aes(cat,ratio,fill=cat))+
  geom_text(aes(label=round(ratio,2),vjust=-0.5),size=6)+
  xlab(NULL)+
  ggtitle("Klf4 OE DEgene in LRgene enrichment")+
  geom_bar(stat="identity")+
  geom_signif(comparisons = list(c("Total","LRgene"),
                                 c("Total","eLR"),
                                 c("Total","Klf4_target_eLR")),
              annotations = as.character(p2star(pvalue.list)),
              y_position = c(0.55,0.65,0.75),
              size = 1,
              textsize = 6)+
  scale_fill_manual(values = c("grey",blues(4)[-1]))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.8),breaks = seq(0,0.8,0.1))+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")
ggsave(filename = myFileName(prefix = "res/fig/Fig6_Klf4_KO_DEgene_in_LRgene",suffix = ".png"),
       width = 8,height = 8,dpi = 350,units = "in",bg = "white")


#####-------3.4 venndigram-------

tmp.list <- list(DE = DE_gene,
                 LR = LRgene.1,
                 eLR = eLR,
                 Klf4_target_eLR = Klf4_eLR_target)
png(filename = myFileName(prefix = "res/fig/Fig6_Klf4_OE_DE_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = T,
     label=T,
     fill=c("white",blues(4)[-1]),
     col="black",
     bg="grey")
dev.off()

#####-----3.5 barplot revised --------


pvalue.list <- NULL
pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(DE_gene,mm10_gene.set)),
                                                  N = length(mm10_gene.set),
                                                  k = length(intersect(LRgene.1,DE_gene)),
                                                  n = length(LRgene.1),
                                                  lower.tail = F))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.1,DE_gene)),
                                                  N = length(LRgene.1),
                                                  k = length(intersect(eLR,DE_gene)),
                                                  n = length(eLR),
                                                  lower.tail = T))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(eLR,DE_gene)),
                                                  N = length(eLR),
                                                  k = length(intersect(Klf4_eLR_target,DE_gene)),
                                                  n = length(Klf4_eLR_target),
                                                  lower.tail = T))


pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(DE_gene,mm10_gene.set)),
                                                  N = length(mm10_gene.set),
                                                  k = length(intersect(eLR,DE_gene)),
                                                  n = length(eLR),
                                                  lower.tail = F))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(DE_gene,mm10_gene.set)),
                                                  N = length(mm10_gene.set),
                                                  k = length(intersect(Klf4_eLR_target,DE_gene)),
                                                  n = length(Klf4_eLR_target),
                                                  lower.tail = F))

ggplot(data = data.plot,aes(cat,ratio,fill=cat))+
  geom_text(aes(label=round(ratio,2),vjust=-0.5),size=12)+
  xlab(NULL)+
  ggtitle("Klf4 OE DEgene \n in  LRgene enrichment")+
  geom_bar(stat="identity")+
  geom_signif(comparisons = list(c("Total","LRgene"),
                                 c("LRgene","eLR"),
                                 c("eLR","Klf4_target_eLR"),
                                 c("Total","eLR"),
                                 c("Total","Klf4_target_eLR")),
              annotations = paste0("p=",signif(pvalue.list,3)),
              #y_position = seq(0.55,0.8,length.out=5),
              y_position = c(0.60,0.54,0.49,0.70,0.76),
              size = 2,
              textsize = 10)+
  scale_fill_manual(values = c("grey",blues(4)[-1]))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.9),breaks = seq(0,0.9,0.1))+
  scale_x_discrete(labels=c("all","LR","eLR","KLF4 \n target eLR"))+
  theme_cowplot(font_size = 30)+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")
ggsave(filename = myFileName(prefix = "res/fig/Fig6_Klf4_OE_DEgene_in_LRgene_value",
                             suffix = ".png"),
       width = 8,height = 8,dpi = 350,units = "in",bg = "white")


###------ 4. kLF binding sites----------

###-------4.1 load data--------

### mm9 gene

mm9KG_txdb <- makeTxDbFromGFF(file = "D://Ubuntu/wlt/igenomes/Mus_musculus/UCSC/mm9/Annotation/genes.gtf")
mm9.gene <- genes(mm9KG_txdb)
mm9.gene.list <- mm9.gene$gene_id
mm9.tss <- promoters(mm9.gene,upstream = 2000, downstream = 2001)
mm9.tss


### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)



### TFBS
tmp.dir <- "res/meme/res/TFbindsites/"
#list.files(path = tmp.dir)
tmp.samples <- c("2cell_early","2cell","4cell","8cell","icm")
tmp.files.1 <- paste0(tmp.dir,tmp.samples,"_peaks_klf4_bindsite.bed")
tmp.files.2 <- paste0("res/meme/res/bed/",tmp.samples,"_peaks.bed")

tmp.res.list <- lapply(1:length(tmp.samples),FUN = function(ii){
  tmp.klf4.sites <- import.bed(tmp.files.1[ii])  
  tmp.ATAC.sites <- import.bed(tmp.files.2[ii])
  ### TFBS
  Klf4_biding_gene <- subsetByOverlaps(mm9.tss,tmp.klf4.sites,type = "any")
  ATAC_gene <- subsetByOverlaps(mm9.tss,tmp.ATAC.sites,type = "any")
  
  #export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
  M <- length(Klf4_biding_gene)
  N <- length(ATAC_gene)
  
  tmp.ATAC.gene <- subsetByOverlaps(mm9.tss[LRgene,],tmp.ATAC.sites,type = "any")
  tmp.Klf4_ATAC.gene <- subsetByOverlaps(tmp.ATAC.gene,Klf4_biding_gene)
  
  
  k <- length(tmp.Klf4_ATAC.gene)
  n <- length(tmp.ATAC.gene)
  
  x <- M/N
  y <- k/n
  p <- phyper(k-1,M,N-M,n,lower.tail = F)
  #tmp.df <- data.frame(expected = x, obs = y,pvalue = p,sample=tmp.samples[ii],stringsAsFactors = F)
  tmp.df <- data.frame(sample=tmp.samples[ii],
                       value=c(x,y),
                       pvalue=c(p,p),
                       group=c("expected","observed"),
                       set = c("all","LRgene"),
                       stringsAsFactors = F)
  return(tmp.df)
})


tmp.res.df <- Reduce(rbind,tmp.res.list)

tmp.res.df$sample <- factor(tmp.res.df$sample,
                            levels = c("2cell_early","2cell","4cell","8cell","icm"))
tmp.res.df.1 <- tmp.res.df

####------4.2 eLR---------
tmp.df <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
LRpairs <- tmp.df %>%
  dplyr::filter(abs(PCC) > 0.1) %>%
  pull(LRpairs)

Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)

#eLR <- LRgene
### TFBS
tmp.dir <- "res/meme/res/TFbindsites/"
#list.files(path = tmp.dir)
tmp.samples <- c("2cell_early","2cell","4cell","8cell","icm")
tmp.files.1 <- paste0(tmp.dir,tmp.samples,"_peaks_klf4_bindsite.bed")
tmp.files.2 <- paste0("res/meme/res/bed/",tmp.samples,"_peaks.bed")

tmp.res.list <- lapply(1:length(tmp.samples),FUN = function(ii){
  tmp.klf4.sites <- import.bed(tmp.files.1[ii])  
  tmp.ATAC.sites <- import.bed(tmp.files.2[ii])
  ### TFBS
  Klf4_biding_gene <- subsetByOverlaps(mm9.tss,tmp.klf4.sites,type = "any")
  ATAC_gene <- subsetByOverlaps(mm9.tss,tmp.ATAC.sites,type = "any")
  
  #export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
  M <- length(Klf4_biding_gene)
  N <- length(ATAC_gene)
  
  tmp.ATAC.gene <- subsetByOverlaps(mm9.tss[LRgene,],tmp.ATAC.sites,type = "any")
  tmp.Klf4_ATAC.gene <- subsetByOverlaps(tmp.ATAC.gene,Klf4_biding_gene)
  
  
  k <- length(tmp.Klf4_ATAC.gene)
  n <- length(tmp.ATAC.gene)
  
  x <- M/N
  y <- k/n
  p <- phyper(k-1,M,N-M,n,lower.tail = F)
  #tmp.df <- data.frame(expected = x, obs = y,pvalue = p,sample=tmp.samples[ii],stringsAsFactors = F)
  tmp.df <- data.frame(sample=tmp.samples[ii],
                       value=y,
                       pvalue=p,
                       group="observed",
                       set = "eLRgene",
                       stringsAsFactors = F)
  return(tmp.df)
})


tmp.res.df <- Reduce(rbind,tmp.res.list)

tmp.res.df$sample <- factor(tmp.res.df$sample,
                            levels = c("2cell_early","2cell","4cell","8cell","icm"))

tmp.res.df.2 <- tmp.res.df



####------4.3 Klf4 target gene---------
tmp.df <- read.delim(file = "res/txt/eLRgene_network_20210512.txt",stringsAsFactors = F,header = T)
#head(tmp.df)
tmp.gene <- tmp.df %>%
  dplyr::filter(source == "Klf4") %>%
  pull(target) %>%
  unique()
LRgene <- intersect(tmp.gene,mm9.gene.list)

#eLR <- LRgene
### TFBS
tmp.dir <- "res/meme/res/TFbindsites/"
#list.files(path = tmp.dir)
tmp.samples <- c("2cell_early","2cell","4cell","8cell","icm")
tmp.files.1 <- paste0(tmp.dir,tmp.samples,"_peaks_klf4_bindsite.bed")
tmp.files.2 <- paste0("res/meme/res/bed/",tmp.samples,"_peaks.bed")

tmp.res.list <- lapply(1:length(tmp.samples),FUN = function(ii){
  tmp.klf4.sites <- import.bed(tmp.files.1[ii])  
  tmp.ATAC.sites <- import.bed(tmp.files.2[ii])
  ### TFBS
  Klf4_biding_gene <- subsetByOverlaps(mm9.tss,tmp.klf4.sites,type = "any")
  ATAC_gene <- subsetByOverlaps(mm9.tss,tmp.ATAC.sites,type = "any")
  
  #export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
  M <- length(Klf4_biding_gene)
  N <- length(ATAC_gene)
  
  tmp.ATAC.gene <- subsetByOverlaps(mm9.tss[LRgene,],tmp.ATAC.sites,type = "any")
  tmp.Klf4_ATAC.gene <- subsetByOverlaps(tmp.ATAC.gene,Klf4_biding_gene)
  
  
  k <- length(tmp.Klf4_ATAC.gene)
  n <- length(tmp.ATAC.gene)
  
  x <- M/N
  y <- k/n
  p <- phyper(k-1,M,N-M,n,lower.tail = F)
  #tmp.df <- data.frame(expected = x, obs = y,pvalue = p,sample=tmp.samples[ii],stringsAsFactors = F)
  tmp.df <- data.frame(sample=tmp.samples[ii],
                       value=y,
                       pvalue=p,
                       group="observed",
                       set = "Klf4eLRgene",
                       stringsAsFactors = F)
  return(tmp.df)
})


tmp.res.df <- Reduce(rbind,tmp.res.list)

tmp.res.df$sample <- factor(tmp.res.df$sample,
                            levels = c("2cell_early","2cell","4cell","8cell","icm"))

tmp.res.df.3 <- tmp.res.df


tmp.res.df <- rbind(tmp.res.df.1,tmp.res.df.2,tmp.res.df.3)

tmp.res.df$set <- factor(tmp.res.df$set,levels = c("all","LRgene","eLRgene","Klf4eLRgene"))
p <- ggplot(data = tmp.res.df,aes(x=sample,y=value,color=set))+
  geom_path(aes(group=set),size=1.5)+
  geom_point(size=6)+
  scale_color_manual(values = c("darkgrey",blues(9)[c(6,8,9)]))+
  scale_x_discrete(labels=c("2cell_early","2cell","4cell","8cell","ICM"))+
  scale_y_continuous(limits=c(0.4,1),breaks = seq(0.4,1,0.1),labels =  seq(0.4,1,0.1))+
  labs(x = NULL,y="ratio",color=NULL)+
  theme_cowplot(font_size = 28)+
  theme(axis.line = element_line(size = 2),
        axis.text.x = element_text(size = 28,angle = 45,vjust = 1,hjust = 1),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(size = 1.5),legend.position = c(0.1,0.8))
p
ggsave(plot = p,filename = myFileName(prefix = "res/fig/fig6_Klf4_binding_enrichment",suffix = ".png"),
       width = 8,height = 8,dpi = 350)

#####-----------5.Rcistarget---------------

####---5.1 basic use------

data("motifAnnotations_mgi")
motifRankings <- importRankings("E://wlt/project/data/cistarget/mm9-tss-centered-10kb-10species.mc9nr.feather")

eLRgene <- read.delim(file = "res/homer/20210511/mm9.tss.eLRgene_sorted.bed",stringsAsFactors = F,header = F)

geneLists <- list(eLRgene=eLRgene$V4)
motifEnrichmentTable_wGenes <- cisTarget(geneLists, 
                                         motifRankings = motifRankings,
                                         motifAnnot = motifAnnotations_mgi,
                                         nesThreshold = 0, 
                                         geneErnMethod="aprox")

tmp.weblogo <- addLogo(motifEnrichmentTable_wGenes)
tmp.weblogo$logo[1]
data.plot <- tmp.weblogo[1:5,]
data.plot$TF_highConf <- gsub(pattern = " .*",replacement = "",x = data.plot$TF_highConf)
data.plot$TF_highConf <- factor(data.plot$TF_highConf,levels = rev(data.plot$TF_highConf))
ggplot(data.plot,aes(TF_highConf,NES))+
  geom_bar(stat="identity",fill="#75aadb")+
  geom_text(aes(label=NES),size=12,hjust=0)+
  coord_flip()+
  theme_cowplot(font_size = 28)+
  scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = seq(0,5,by=1))+
  scale_x_discrete(expand = c(0,0))+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))

ggsave(filename = myFileName(prefix = "res/fig/Fig6_RcisTarget",suffix = ".png"),
       width = 8,height = 4,dpi = 350)


ggplot(data.plot,aes(TF_highConf,NES))+
  geom_bar(stat="identity",fill="#75aadb")+
  geom_text(aes(label=NES),size=12,hjust=0)+
  coord_flip()+
  ggprism::theme_prism(base_size = 28,base_family = "serif",base_line_size = 2)+
  scale_y_continuous(expand = c(0,0),limits = c(0,5),breaks = seq(0,5,by=1))+
  scale_x_discrete(expand = c(0,0))

?ggprism::theme_prism

library(DT)
datatable(tmp.weblogo[1:5,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE)


#####--------5.2 plot logo-----------

library(universalmotif)
library(ggseqlogo)
data(ggseqlogo_sample)
pfms_dna$MA0018.2

test <- read_meme(file = "database/motif/Klf6_jaspar.meme")
view_motifs(test)

test <- read.table(file = "database/motif/GPD1.output.txt",stringsAsFactors = F)

tmp.file <- "database/motif/GPD1.output.txt"
tmp_load_motif <- function(tmp.file,tmp.name="Gdp1"){
  tmp.data <- read.table(file = tmp.file,stringsAsFactors = F)
  #rownames(tmp.data) <- tmp.data[,1]
  tmp.data <- tmp.data[,-1]
  tmp.motif <- create_motif(input = as.matrix(tmp.data),
               type = "PCM",
               alphabet = "DNA",
               name = tmp.name)
  #tmp.motif <- convert_type(tmp.motif,type = "PPM")
  return(tmp.motif)
}

Klf4.motif <- read_meme(file = "database/motif/Klf4_jaspar.meme")
Gpd1.motif <- tmp_load_motif(tmp.file = "database/motif/GPD1.output.txt",
                             tmp.name = "Gpd1")

DNA.motif <- create_motif(mat, alphabet = "DNA")

Spats2.motif <- tmp_load_motif(tmp.file = "database/motif/SPATS2.output.txt",
                               tmp.name = "Spats2")
Trove2.motif <- tmp_load_motif(tmp.file = "database/motif/TROVE2.output.txt",
                               tmp.name = "Trove2")
Klf6.motif <- read_meme(file = "database/motif/Klf6_jaspar.meme")

tmp.plot.list <- list()
p <- view_motifs(Klf4.motif) +
  theme_map()
tmp.plot.list <- c(tmp.plot.list,list(p))
p <- view_motifs(Gpd1.motif) +
  theme_map()
tmp.plot.list <- c(tmp.plot.list,list(p))
p <- view_motifs(Spats2.motif) +
  theme_map()
tmp.plot.list <- c(tmp.plot.list,list(p))
p <- view_motifs(Trove2.motif) +
  theme_map()
tmp.plot.list <- c(tmp.plot.list,list(p))
p <- view_motifs(Klf6.motif) +
  theme_map()
tmp.plot.list <- c(tmp.plot.list,list(p))

cowplot::plot_grid(plotlist = tmp.plot.list,
                   ncol = 1,align = "hv")


library(patchwork)
view_motifs(motifs = c(Klf4.motif,
                       Gpd1.motif,
                       Spats2.motif,
                       Trove2.motif,
                       Klf6.motif),
            min.overlap = 10,
            tryRC = T) & 
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank())
ggsave(filename = myFileName(prefix = "res/fig/fig6_motif",
                             suffix = ".pdf"),
       width = 8,height = 6,dpi = 350)




