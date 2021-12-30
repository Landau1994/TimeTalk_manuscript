####Fig3
####20211003
####count Ligand receptor gene ratio
####edit this in 20211021
rm(list = ls())
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
library(ggrepel)
library(VennDiagram)
library(nichenetr)
library(igraph)
library(DiagrammeR)
library(DiagrammeRsvg)
library(patchwork)
library(GenomicFeatures)
source(file = "code/myUtils.R")
library(readxl)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


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
strike_freedom_gundam_color <- c("#f6383e","#eff3ff","#1f55d9","#373b35","#b1b0c2","#f2be58")
scales::show_col(strike_freedom_gundam_color)


scales::show_col(blues(10))
scales::show_col(ylord(10))
###hic.pca.red <- colorRampPalette(c("#1a469c","black","#e7141a"))
hic.pca.red <- colorRampPalette(c("blue","gray1","red"))
hic.pca.redwhite <- colorRampPalette(c("#1d1856","navyblue","white","red4","#861617"))
hic.pca.orange <- colorRampPalette(c("#2f2583","black","#f9b232"))
hic.pca.skyblue <- colorRampPalette(c("skyblue","black","orange"))
hic.pca.gundam <- colorRampPalette(strike_freedom_gundam_color[1:3])

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
divergentcolor.gundam <- function(n){
  if(n <= length(strike_freedom_gundam_color)){
    colors <- strike_freedom_gundam_color[1:n]
  }
  else{
    colors <- (grDevices::colorRampPalette(strike_freedom_gundam_color))(n)
  }
}

navy <- colorRampPalette(colors = c("#3e56a6", "#fffbfb", "#ee252a"))
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
coldwarm <- colorRampPalette(colors = mypalette)

bkcolor <- c(colorRampPalette(c(brewer.pal(9,"Blues" )[4:9],"#1a1919"))(50),
             colorRampPalette(c("#1a1919",rev(brewer.pal( 9,"YlOrBr" ))[1:6]))(50))

bkcolor <- colorRampPalette(colors = bkcolor)
mycolor.bar(bkcolor(100),min = -1)
mycolor.bar(blues(10),min = -1)
coffeshop <- colorRampPalette(colors =c("#dd9b53","#d25b12","seagreen3","#c2a38a","#b8a896","grey"))
scales::show_col(colours = coffeshop(8))
scales::show_col(divergentcolor(15))


scales::show_col(colorRampPalette(c("#00bcd4","lightgrey","#ffeb3b"))(5))


shaman.color <- c(c("#08306b","#084d96","#3f51b5","#2196f3","#03a9f4","#00bcd4"),
                  colorRampPalette(c("#00bcd4","lightgrey","#ffeb3b"))(5)[2:4],
                  rev(c("#800026","#B60026","#ff5722","#ff9800","#ffc107","#ffeb3b")))

shaman.color <- colorRampPalette(colors = shaman.color)

par(mfrow=c(2,2))
mycolor.bar(shaman.color(100),min = -1,vertical = F)
mycolor.bar(rdbu(100),min = -1,vertical = F)
mycolor.bar(coldwarm(100),min = -1,vertical = F)
mycolor.bar(solarExtra(100),min = -1,vertical = F)
dev.off()
scales::show_col(solarExtra(5))
scales::show_col(shaman.color(11))


tmp.color <- rev(c("#f4ea0c","#ffc807",
                   "#f6901f","#f05e23",
                   "#ed3324","#ec1c25",
                   "#bc243c","#8d3966",
                   "#5a5391","#366cb3",
                   "#1f65af","#174c81",
                   "#163256","#0e192a",
                   "#000000"))
bkro <- colorRampPalette(colors = tmp.color)
mycolor.bar(hic.pca.gundam(100),min = -1)

####-----------1. visualization of eLR screen---------------------


tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")

####plot
cutoff <- 0.1
data.plot <- tmp.res.cor.1 %>%
  mutate(Rank=1:nrow(tmp.res.cor.1)) %>%
  mutate(group=ifelse(PCC > cutoff,"p-eLR",ifelse(PCC < -cutoff,"n-eLR","not sig"))) %>%
  mutate(group=factor(group,levels = c("p-eLR","not sig","n-eLR"))) %>%
  mutate(LRpairs=gsub(x = LRpairs,pattern = "_",replacement = "-"))




tmp.show <- 5
tmp.label.1 <- data.plot %>%
  head(tmp.show) %>%
  pull(LRpairs)

tmp.label.2 <- data.plot %>%
  tail(tmp.show) %>%
  pull(LRpairs)

tmp.label <- c(tmp.label.1,tmp.label.2)


p <- ggplot(data.plot,aes(Rank,PCC,color=group,label=ifelse(LRpairs %in% tmp.label,LRpairs,"")))+
  geom_point(alpha=0.8)+
  geom_text_repel(box.padding = 0.5,
                  max.overlaps = Inf)+
  scale_color_manual(values = c("#c00000","grey","#0070c0"))+
  guides(color = guide_legend(override.aes = aes(label = "")))+
  scale_x_continuous(breaks = seq(0,1200,100),labels = seq(0,1200,100))+
  scale_y_continuous(breaks = seq(-1,1,0.2),labels = seq(-1,1,0.2))+
  geom_hline(yintercept = cutoff,linetype="dotdash")+
  geom_hline(yintercept = -cutoff,linetype="dotdash")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1))
p
ggsave(p,filename = myFileName(prefix = "res/fig/Fig3_eLR_show_top5",suffix = ".jpg"),
       width = 9,height = 9,dpi = 350)


tmp.label <- c("Fgf4-Fgfr1","Fgf4-Fgfr2","Bmp4-Bmpr2",
               "Gdf9-Bmpr1b","Ntn1-Dcc","Ntn1-Unc5b",
               "Rspo2-Lgr5","Rspo2-Lgr6",'Hbegf-Erbb4',
               "Pdgfc-Pdgfra","Igf2-Igf2r",
               "Ctsd-Lrp1")

p <- ggplot(data.plot,aes(Rank,PCC,color=group,label=ifelse(LRpairs %in% tmp.label,LRpairs,"")))+
  geom_point(alpha=0.8,size=6)+
  geom_text_repel(box.padding = 0.5,
                  size=6,fontface = "italic",
                  max.overlaps = Inf)+
  scale_color_manual(values = c("#c00000","grey","#0070c0"),name=NULL)+
  guides(color = guide_legend(override.aes = aes(label = "")))+
  scale_x_continuous(breaks = seq(0,1200,100),labels = seq(0,1200,100))+
  scale_y_continuous(breaks = seq(-1,1,0.2),labels = seq(-1,1,0.2),limits = c(-1,1))+
  geom_hline(yintercept = cutoff,linetype="dotdash")+
  geom_hline(yintercept = -cutoff,linetype="dotdash")+
  theme_cowplot(font_size = 30)+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(angle = 45,vjust = 0.5),
        legend.position = "top")
p
ggsave(p,filename = myFileName(prefix = "res/fig/Fig3_eLR_show_known",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


####plot for model
p <- ggplot(data.plot,aes(Rank,PCC,color=group,label=ifelse(LRpairs %in% tmp.label,LRpairs,"")))+
  geom_point(alpha=0.8,size=3)+
  # geom_text_repel(box.padding = 0.5,
  #                 max.overlaps = Inf)+
  scale_color_manual(name=NULL,
                     values = c("#c00000","grey","#0070c0"),
                     labels = c("p-eLR(322)","not sig","n-eLR(132)"))+
  #guides(color = guide_legend(override.aes = aes(label = "")))+
  scale_x_continuous(breaks = seq(0,1200,100),
                     labels = seq(0,1200,100))+
  scale_y_continuous(breaks = seq(-1,1,0.2),
                     labels = seq(-1,1,0.2),
                     limits = c(-1,1))+
  geom_hline(yintercept = cutoff,linetype="dotdash",)+
  geom_hline(yintercept = -cutoff,linetype="dotdash")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1))
p
ggsave(p,filename = myFileName(prefix = "res/fig/Fig3_eLR_show_formodel",
                               suffix = ".pdf"),
       width = 8,height = 8,dpi = 350)





####----2. overlap with ZGA gene------

####----2.1 load data-------


### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
LRgene.all <- LRgene


####eLR gene
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
cutoff <- 0.1
tmp.LR.pairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)
eLR_gene <- c(tmp.L.gene,tmp.R.gene)






# reprogramming_factor_candidate <- read.delim(file = "database/genes_365_for_iPSCs_202001303_DBTMEE.tsv",
#                                              sep = "\t",stringsAsFactors = F) %>%
#   pull("Gene")
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)
genelist <- list(eLR_gene=eLR_gene,
                 zga_gene=zga_gene_zhangyi)

####-----2.2 plot-----------

###draw venn diagram
scales::show_col(colours = shaman.color(20))
scales::show_col(colours = c("#84c683","#fec04e"))
scales::show_col(colours = c("#E69F00","#009E73"))
tmp.color <- c("#80bda0","#e0ac61")
scales::show_col(ggthemes::colorblind_pal()(9))


genelist <- list(zga_gene=zga_gene_zhangyi,
                eLR_gene=eLR_gene)

#### get pvalue
library(EnsDb.Mmusculus.v79)
mm10_gene <- genes(EnsDb.Mmusculus.v79)
mm10_gene <- as.data.frame(mm10_gene)
mm10_gene <- mm10_gene %>%
  dplyr::filter(gene_biotype == "protein_coding")
mm10_gene.set <- mm10_gene$symbol

k <- length(intersect(eLR_gene,zga_gene_zhangyi))
M <- length(zga_gene_zhangyi)
N <- length(mm10_gene.set)
n <- length(eLR_gene)
x <- M/N
y <- k/n
pvalue <- phyper(k-1,M,N-M,n,lower.tail = F)

tmp.df <- data.frame(value=c(x,y),group=c("expected","observed"))
ggplot(tmp.df,aes(group,value,fill=group))+
  geom_bar(stat = "identity")+
  ggsignif::geom_signif(comparisons = list(c("expected","observed")),
                        annotations = paste0("p=",round(pvalue,6)),
                        tip_length = 0.3,
                        size = 2,y_position = 0.15,
                        textsize = 10)+
  xlab(NULL)+
  ylab("ratio")+
  scale_fill_manual(values = c("#80bda0","#e0ac61"))+
  scale_x_discrete(labels=c("ZGA gene","eLR gene"))+
  scale_y_continuous(limits = c(0,0.18),expand = c(0,0))+
  theme(legend.position = "none")+
  theme_cowplot(font_size = 45)+
  theme(axis.line = element_line(size = 2),
        axis.text.x = element_text(angle = 30,hjust = 0.5,vjust = 0.8),
        axis.ticks = element_line(size = 2),
        legend.position = "none")  
ggsave(filename = myFileName(prefix = "res/fig/fig3_bar",suffix = ".png"),width = 6,height = 8,dpi = 350)

png(myFileName(prefix = paste0("res/fig/fig3_venn_zga_gene"),
               suffix = ".png"),units = "in",res = 350,
    width = 8,height = 8)
plot(eulerr::euler(genelist,shape = "ellipse"),
     quantities = list(fontsize=30),
     label=list(fontsize=30),
     fill=tmp.color,
     main=list(label=paste0("p=",sprintf("%0.3f",pvalue)),fontsize=24),
     col="black")
dev.off()
?euler

pdf(myFileName(prefix = paste0("res/fig/fig3_venn_zga_gene"),
               suffix = ".pdf"),
    width = 8,height = 8)
plot(eulerr::euler(genelist,shape = "ellipse"),
     quantities = list(fontsize=30),
     label=list(fontsize=30),
     fill=tmp.color,
     main=list(label=paste0("p=",sprintf("%0.3f",pvalue)),fontsize=24),
     col="black")
dev.off()




#####--------3. nichenet----------



######---------3.2.1 load data------------
early.scRNAseq.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
early.scRNAseq.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")

ligand_target_matrix = readRDS("res/R/nichenet_ligand_target_matrix.rds")
lr_network <- readRDS(file = "res/R/nichenet_lr_network_mouse.rds")
mouse_sig_network <- readRDS(file = "res/R/nichenet_mouse_sig_network.rds")
weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")
ligand_tf_matrix = readRDS(file = "res/R/nichenet_ligand_tf_matrix.rds")

#####---------3.2.2 define expressed genes---------
tmp.id <- early.scRNAseq.meta %>%
  dplyr::filter(Stage == "early2cell") %>%
  pull(cell_id)
expressed_genes_sender <- early.scRNAseq.exp[,tmp.id] %>%
  apply(2,function(x){2^x-1}) %>% 
  apply(1,function(x){log2( mean(x) +1)}) %>% 
  .[. >= 4] %>% 
  names()
expressed_genes_receiver <- expressed_genes_sender



######-------3.2.3 define gene sets, ligands------
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene) 
geneset_oi <- zga_gene_zhangyi
background_expressed_genes = expressed_genes_receiver  %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% 
  dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()







########--------3.2.4 Perform NicheNet’s ligand activity analysis on the gene set of interest--------
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities %>% arrange(-pearson) 
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
best_upstream_ligands
# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity


######------------3.2.5 Infer target genes of top-ranked ligands and visualize in a heatmap ----------

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0)
nrow(active_ligand_target_links)
## [1] 143
head(active_ligand_target_links)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands",
                      "ZGA gene", 
                      color = "purple",
                      x_axis = T,
                      legend_position = "top", 
                      x_axis_position = "top",
                      legend_title = "Regulatory potential") + 
  scale_fill_gradientn(colours = hic.orage(256),
                        breaks = c(0,0.009),
                        labels = c(0,0.009))+ 
  theme(text = element_text(size = 24),
        axis.text.x.top= element_blank(),
        axis.ticks.x = element_blank())

p_ligand_target_network
round(range(vis_ligand_target),2)
ggsave(filename = myFileName(prefix = "res/fig/Fig3_ZGA_ligand_nichenet",suffix = ".png"),
       width = 24,height = 8,dpi = 350)



#####-----------3.2.6 : Ligand-receptor network inference for top-ranked ligands------------------------


# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% 
  dplyr::filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% 
  dplyr::filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% 
  dplyr::select(-to) %>% 
  as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized CAF-ligands","Receptors expressed by malignant cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network


ligand_pearson_matrix = ligand_activities %>% 
  dplyr::select(pearson) %>% 
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", 
                                                              x_axis_position = "top", legend_title = "activity")
p_ligand_pearson = p_ligand_pearson+scale_fill_gradientn(colors  = blues(256),
                                      breaks = round(range(vis_ligand_pearson),2),
                                      labels = round(range(vis_ligand_pearson),2))

#?make_heatmap_ggplot

# (p_ligand_pearson | p_ligand_target_network) / p_ligand_receptor_network
# ggsave(filename = "res/fig/nichenet_res_test_early2cell_20210130.jpg",width = 8,height = 8,dpi = 350)

plot_grid(p_ligand_pearson+ylab("prioritized Ligands")+theme(text = element_text(size = 30),
                                                             legend.text = element_text(size = 24),
                                                             legend.text.align = 0.5,
                                                             axis.text.x.top = element_blank(),
                                                             axis.ticks.x = element_blank()),
          p_ligand_target_network+ylab(NULL)+
            theme(axis.text.y = element_blank(),
                  legend.text = element_text(size = 24),
                  axis.ticks.y = element_blank()),rel_widths = c(2,8))
ggsave(filename = myFileName(prefix = "res/fig/fig3_nichenet_res_ZGA",suffix = ".png"),
       width = 16,height = 9,dpi = 350,units = "in")

tmp.mat <- vis_ligand_pearson
ComplexHeatmap::pheatmap(mat = tmp.mat[rev(1:nrow(tmp.mat)),],
                         color = rdbu(100),
                         cluster_rows = F,cluster_cols = F)


######----------3.2.7 get graph------
ligands_all = c("Tgfb2") # this can be a list of multiple ligands if required
targets_all = c("Cdkn1a")

ligand_target_matrix = readRDS("res/R/nichenet_ligand_target_matrix.rds")
lr_network <- readRDS(file = "res/R/nichenet_lr_network_mouse.rds")
mouse_sig_network <- readRDS(file = "res/R/nichenet_mouse_sig_network.rds")
weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")
ligand_tf_matrix = readRDS(file = "res/R/nichenet_ligand_tf_matrix.rds")
head(lr_network)
# weighted_networks = construct_weighted_networks(lr_network, mouse_sig_network, mouse_gr_network,
#                                                 source_weights_df)

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, 
                                                     ligands_all = ligands_all, 
                                                     targets_all = targets_all, 
                                                     weighted_networks = weighted_networks)
# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 5.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 5.75)
graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network, 
                                                  ligands_all = ligands_all, 
                                                  targets_all = targets_all, 
                                                  sig_color = "indianred", 
                                                  gr_color = "steelblue")
unique(active_signaling_network_min_max$sig)
# To render the graph: uncomment following line of code
test <- DiagrammeR::render_graph(graph_min_max, layout = "tree",output = "graph")
test
export_svg(test) %>%
  charToRaw %>% rsvg %>% ('res/fig/test_nichenet_path_width.png')
test.tttt <- export_svg(test) %>%
  charToRaw


test <- graph_min_max %>% to_igraph() 

ttt.layout = layout_with_sugiyama(test)

tmp.vertext <- vertex.attributes(test)
tmp.edge <- edge.attributes(test)
tmp.edge
tmp.layout <- ttt.layout$layout
jpeg(filename = myFileName(prefix = "res/fig/fig3_test_nichenet_path_width",suffix = ".jpg"),
     width = 8,height = 8,res = 350,units = "in")
plot.igraph(test,
            layout=tmp.layout,
            vertex.size = 30,
            vertex.label.font=18,
            vertex.label.cex=1,
            vertex.frame.color=tmp.vertext$fillcolor,
            vertex.color=tmp.vertext$fillcolor,
            vertex.label.color=tmp.vertext$fontcolor,
            edge.width = 10*tmp.edge$penwidth)
dev.off()


pdf(file  = myFileName(prefix = "res/fig/fig3_test_nichenet_path_width",suffix = ".pdf"),
    width = 8,height = 8)
plot.igraph(test,
            layout=tmp.layout,
            vertex.size = 30,
            vertex.label.font=18,
            vertex.label.cex=1,
            vertex.frame.color=tmp.vertext$fillcolor,
            vertex.color=tmp.vertext$fillcolor,
            vertex.label.color=tmp.vertext$fontcolor,
            edge.width = 10*tmp.edge$penwidth)
dev.off()


###-------4.gene sets overlap---------
mm9KG_txdb <- makeTxDbFromGFF(file = "D://Ubuntu/wlt/igenomes/Mus_musculus/UCSC/mm9/Annotation/genes.gtf")
mm9.gene <- genes(mm9KG_txdb)
mm9.gene.list <- mm9.gene$gene_id

test <- transcripts(mm9KG_txdb)
test <- exons(mm9KG_txdb)
??GenomicFeatures

columns(mm9KG_txdb)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
columns(TxDb.Mmusculus.UCSC.mm9.knownGene)



seqinfo(mm9KG_txdb)
seqlengths(mm9KG_txdb)



###-------4.1 essential gene---------

###-------4.1.1 LR--------

### load essential gene list
tmp.dir <- "database/deg_annotation_e.csv/deg_annotation_e.csv"
tmp.df <- read.delim(file = tmp.dir,sep = ";",stringsAsFactors = F,header = F)
essential.gene <- tmp.df %>%
  dplyr::filter(V8 == "Mus musculus") %>%
  pull(V3)
head(essential.gene)
essential.gene <- intersect(essential.gene,mm9.gene.list)

### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
LRgene.all <- LRgene

#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(essential.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,essential.gene))
n <- length(LRgene)
x <- M/N
y <- k/n
p <- phyper(k-1,M,N-M,n,lower.tail = F)
tmp.df.1 <- data.frame(value = c(x,y),group=c("all","LR"),stringsAsFactors = F)
pvalue.list <- p

####-------4.1.2 eLR------------------
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
eLRgene <- LRgene
#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(essential.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,essential.gene))
n <- length(LRgene)
x <- M/N
y <- k/n


pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.all,essential.gene)),
                                                  N = length(LRgene.all),
                                                  k = length(intersect(eLRgene,essential.gene)),
                                                  n = length(eLRgene)))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(mm9.gene.list,essential.gene)),
                                                  N = length(mm9.gene.list),
                                                  k = length(intersect(eLRgene,essential.gene)),
                                                  n = length(eLRgene)))


tmp.df.2 <- data.frame(value = y,
                       group=c("eLR"),
                       stringsAsFactors = F)
tmp.df <- rbind(tmp.df.1,tmp.df.2)
tmp.df$group <- factor(tmp.df$group,levels = c("all","LR","eLR"))

ggplot(tmp.df,aes(group,value,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(value,2),vjust=-0.5),size=12)+
  geom_signif(comparisons = list(c("all","LR"),
                                 c("LR","eLR"),
                                 c("all","eLR")),
              annotations = paste0("p=",signif(pvalue.list,4)),
              y_position = c(0.4,0.55,0.65),size = 2,textsize = 10,vjust = - 0.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.75),breaks = seq(0,0.7,by=0.1))+
  labs(x = NULL,y="ratio")+
  ggtitle("Essential gene enrichment")+
  scale_fill_manual(values = c("darkgrey",blues(3)[2:3]))+
  theme_cowplot(font_size = 40)+
  theme(legend.position = "none",
        plot.title = element_text(size = 36,hjust = 0.5),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
ggsave(filename = myFileName(prefix = "res/fig/Fig3_essential_gene_enrichment",
                             suffix = ".png"),width = 8,height = 8,dpi = 350)

#####---------4.1.3 euler plot---------
tmp.list <- list(all = essential.gene,
                 LR = LRgene.all,
                 eLR = eLRgene)
png(filename = myFileName(prefix = "res/fig/Fig3_essential_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Essential gene overlap",fontsize = 24),
     bg="grey")
dev.off()

pdf(file = myFileName(prefix = "res/fig/Fig3_essential_gene_euler",suffix = ".pdf"),
    width = 8,height = 8,bg = "grey")
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Essential gene overlap",fontsize = 24),
     bg="grey")
dev.off()

###-------4.2 housekeeping gene---------

###-------4.2.1 LR--------

### load HK gene
tmp.dir <- "database/HRT_housekeeping_gene_database/Housekeeping_GenesMouse.csv"
tmp.df <- read.delim(file = tmp.dir,sep = ";",stringsAsFactors = F,header = T)
HK.gene <- tmp.df %>%
  pull(Gene) %>%
  unique()
HK.gene <- intersect(HK.gene,mm9.gene.list)

### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
LRgene.all <- LRgene

#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(HK.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,HK.gene))
n <- length(LRgene)
x <- M/N
y <- k/n
p <- phyper(k-1,M,N-M,n,lower.tail = T)
tmp.df.1 <- data.frame(value = c(x,y),group=c("all","LR"),stringsAsFactors = F)
pvalue.list <- p

####-------4.2.2 eLR------------------
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
eLRgene <- LRgene
#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(HK.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,HK.gene))
n <- length(LRgene)
x <- M/N
y <- k/n

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.all,HK.gene)),
                                                  N = length(LRgene.all),
                                                  k = length(intersect(eLRgene,HK.gene)),
                                                  n = length(eLRgene)))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(mm9.gene.list,HK.gene)),
                                                  N = length(mm9.gene.list),
                                                  k = length(intersect(eLRgene,HK.gene)),
                                                  n = length(eLRgene),lower.tail = T))


tmp.df.2 <- data.frame(value = y,
                       group=c("eLR"),
                       stringsAsFactors = F)
tmp.df <- rbind(tmp.df.1,tmp.df.2)
tmp.df$group <- factor(tmp.df$group,levels = c("all","LR","eLR"))

ggplot(tmp.df,aes(group,value,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(value,2),vjust=-0.5),size=12)+
  geom_signif(comparisons = list(c("all","LR"),
                                 c("LR","eLR"),
                                 c("all","eLR")),
              annotations = paste0("p=",signif(pvalue.list,4)),
              y_position = c(0.16,0.05,0.2),
              size = 2,textsize = 10,vjust = -0.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25),breaks = seq(0,0.2,by=.1))+
  labs(x = NULL,y="ratio")+
  ggtitle("Housekeeping gene \n enrichment")+
  scale_fill_manual(values = c("darkgrey",blues(3)[2:3]))+
  theme_cowplot(font_size = 40)+
  theme(legend.position = "none",
        plot.title = element_text(size = 36,hjust = 0.5),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
ggsave(filename = myFileName(prefix = "res/fig/Fig3_HK_gene_enrichment",
                             suffix = ".png"),width = 8,height = 8,dpi = 350)

####--------4.2.3 euler plot---------------
tmp.list <- list(all = HK.gene,
                 LR = LRgene.all,
                 eLR = eLRgene)
png(filename = myFileName(prefix = "res/fig/Fig3_HK_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Housekeeping gene overlap",fontsize = 24),
     bg="grey")
dev.off()

pdf(file = myFileName(prefix = "res/fig/Fig3_HK_gene_euler",suffix = ".pdf"),
    width = 8,height = 8,bg = "grey")
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Housekeeping gene overlap",fontsize = 24),
     bg="grey")
dev.off()






####--------4.3 ZGA gene overlap--------------

###-------4.3.1 LR--------

### load ZGA gene
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)
zga_gene_zhangyi <- intersect(zga_gene_zhangyi,mm9.gene.list)

### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
LRgene.all <- LRgene

#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(zga_gene_zhangyi)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,zga_gene_zhangyi))
n <- length(LRgene)
x <- M/N
y <- k/n
p <- phyper(k-1,M,N-M,n,lower.tail = F)
tmp.df.1 <- data.frame(value = c(x,y),group=c("all","LR"),stringsAsFactors = F)
pvalue.list <- p

####-------4.3.2 eLR------------------
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
eLRgene <- LRgene
#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(zga_gene_zhangyi)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,zga_gene_zhangyi))
n <- length(LRgene)
x <- M/N
y <- k/n

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.all,zga_gene_zhangyi)),
                                                  N = length(LRgene.all),
                                                  k = length(intersect(eLRgene,zga_gene_zhangyi)),
                                                  n = length(eLRgene)))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(mm9.gene.list,zga_gene_zhangyi)),
                                                  N = length(mm9.gene.list),
                                                  k = length(intersect(eLRgene,zga_gene_zhangyi)),
                                                  n = length(eLRgene)))


tmp.df.2 <- data.frame(value = y,
                       group=c("eLR"),
                       stringsAsFactors = F)
tmp.df <- rbind(tmp.df.1,tmp.df.2)
tmp.df$group <- factor(tmp.df$group,levels = c("all","LR","eLR"))

ggplot(tmp.df,aes(group,value,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(value,2),vjust=-0.5),size=12)+
  geom_signif(comparisons = list(c("all","LR"),
                                 c("LR","eLR"),
                                 c("all","eLR")),
              annotations = paste0("p=",signif(pvalue.list,4)),
              y_position = c(0.11,0.15,0.18),
              size = 2,textsize = 10,vjust = -0.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25),breaks = seq(0,0.25,by=0.1))+
  labs(x = NULL,y="ratio")+
  ggtitle("ZGA gene  enrichment")+
  scale_fill_manual(values = c("darkgrey",blues(3)[2:3]))+
  theme_cowplot(font_size = 40)+
  theme(legend.position = "none",
        plot.title = element_text(size = 36,hjust = 0.5),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
ggsave(filename = myFileName(prefix = "res/fig/Fig3_ZGA_gene_enrichment",
                             suffix = ".png"),width = 8,height = 8,dpi = 350)

####--------4.3.3 euler plot---------------
tmp.list <- list(all = zga_gene_zhangyi,
                 LR = LRgene.all,
                 eLR = eLRgene)
png(filename = myFileName(prefix = "res/fig/Fig3_ZGA_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="ZGA gene overlap",fontsize = 24),
     bg="grey")
dev.off()


pdf(file = myFileName(prefix = "res/fig/Fig3_ZGA_gene_euler",suffix = ".pdf"),
    width = 8,height = 8,bg = "grey")
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="ZGA gene overlap",fontsize = 24),
     bg="grey")
dev.off()


# ####-----(drpreacted)4.3 gene age---------
# ####20211021
# ####The following code will remove in the future version
# 
# tmp.path <- "database/essential_gene_chart.xlsx"
# readxl::excel_sheets(tmp.path)
# tmp.df <- readxl::read_xlsx(tmp.path,sheet = "Sheet1",skip = 1)
# 
# colnames(tmp.hom)
# tmp.hom <- read.delim(file = "database/HOM_MouseHumanSequence.rpt",stringsAsFactors = F)
# tmp.hg <- tmp.hom %>%
#   dplyr::filter( Common.Organism.Name == "human" ) %>%
#   dplyr::filter( Symbol %in% tmp.df$gene) %>%
#   dplyr::select( HomoloGene.ID,Symbol)
# tmp.mm <- tmp.hom %>%
#   dplyr::filter( Common.Organism.Name == "mouse, laboratory" ) %>%
#   dplyr::filter( HomoloGene.ID %in% tmp.hg$HomoloGene.ID) %>%
#   dplyr::select( HomoloGene.ID,Symbol)
# tmp.hom <- merge(tmp.hg,tmp.mm,by="HomoloGene.ID") 
# colnames(tmp.hom)[2:3] <- c("human","mouse")
# 
# colnames(tmp.df)
# 
# not_eLR_gene <- setdiff(LRgene.all,eLRgene)
# #### obtain gene_age table
# gene_age.df <- tmp.df %>%
#   dplyr::select(gene,`gene_age(26)`) %>%
#   dplyr::filter(gene %in% tmp.hom$human) %>%
#   mutate(gene=plyr::mapvalues(gene,from = tmp.hom$human,to = tmp.hom$mouse)) %>%
#   mutate(group="not_LR") %>%
#   mutate(group=ifelse(gene %in% eLRgene,"eLR",group)) %>%
#   mutate(group=ifelse(gene %in% not_eLR_gene,"other_LR",group))
# table(gene_age.df$group)
# gene_age.df$`gene_age(26)` <- as.numeric(gene_age.df$`gene_age(26)`)
# ggplot(gene_age.df,aes(group,`gene_age(26)`,fill=group))+
#   geom_boxplot(size=2)+
#   theme_cowplot(font_size = 30)+
#   xlab(NULL)+
#   theme(axis.line = element_line(size = 2),
#         legend.position = "none",
#         axis.ticks = element_line(size = 2))
# ggsave(filename = myFileName(prefix = "res/fig/Fig3_gene_age_LRgene",
#                              suffix = ".png"),width = 8,height = 8,dpi = 350)





####-------5.blastoid analysis------------

####-------5.1 load data-------


####EPS blastoid
tmp.dir <- "data/blastoid/GSE135701/EPS_blastoid/"
tmp.mat <- Read10X(data.dir = tmp.dir)
tttt <- DropletUtils::barcodeRanks(tmp.mat)
br.out <- tttt
tmp.plot <- br.out[!is.na(br.out$fitted),]
tmp.cell.id <- rownames(tmp.plot)
tmp.mat.1 <- tmp.mat[,tmp.cell.id]
seu <- CreateSeuratObject(tmp.mat.1,min.features = 2000,min.cells = 0,project = "EPS_blastoid")
colnames(seu@meta.data)
VlnPlot(seu,features = c("nCount_RNA","nFeature_RNA"))

seu.1 <- seu

####blastocyst
tmp.dir <- "data/blastoid/GSE135701/Blastocyst/"
tmp.mat <- Read10X(data.dir = tmp.dir)
tttt <- DropletUtils::barcodeRanks(tmp.mat)
br.out <- tttt
tmp.plot <- br.out[!is.na(br.out$fitted),]
tmp.cell.id <- rownames(tmp.plot)
tmp.mat.1 <- tmp.mat[,tmp.cell.id]
seu <- CreateSeuratObject(tmp.mat.1,min.features = 2000,min.cells = 0,project = "blastocyst")
colnames(seu@meta.data)
seu.2 <- seu

###merge
seu <- merge(seu.1,seu.2)



seu <- NormalizeData(seu)
seu <- ScaleData(seu)
seu <- FindVariableFeatures(seu,selection.method = "vst",nfeatures = 3000)
seu <- RunPCA(seu,features = VariableFeatures(seu))
# tmp.mat <- as.matrix(seu@assays$RNA@data)
# tmp.genes <- VariableFeatures(seu)
# seu.dist <- 1- cor(tmp.mat[tmp.genes,],method = "pearson")
# seu <- RunTSNE(seu,
#                distance.matrix = seu.dist,
#                dims = 1:20,
#                seed.use = 42,
#                perplixity=100)
seu <- RunTSNE(seu,dims = 1:20,seed.use = 42)
seu <- RunUMAP(seu,dims = 1:20,umap.method = "umap-learn")
seu <- FindNeighbors(seu,dims = c(1:20))
seu <- FindClusters(seu,resolution = 0.2)

UMAPPlot(seu,group.by="orig.ident")

colnames(seu@meta.data)
FeaturePlot(seu,features = c("Ascl2"))



tmp.list <- list(seu.1,seu.2)
names(tmp.list) <- c("EPS_blastoid","blastocyst")
tmp.list <- lapply(names(tmp.list),FUN = function(x){
  tmp.list[[x]] <- NormalizeData(tmp.list[[x]], verbose = FALSE)
  tmp.list[[x]] <- FindVariableFeatures(tmp.list[[x]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
  return(tmp.list[[x]])
})

#reference.list <- pancreas.list[c("PT081", "PT089")]
#?FindIntegrationAnchors
tmp.list <- FindIntegrationAnchors(object.list = tmp.list, 
                                   dims = 1:20,
                                   k.anchor = 5,
                                   k.filter = 30)
#?IntegrateData
seu.integrated <- IntegrateData(anchorset = tmp.list, dims = 1:20)


###RunSeurat pipeline
DefaultAssay(seu.integrated) <- "integrated"
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, npcs = 30, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, reduction = "pca", dims = 1:30,umap.method = "umap-learn")
seu.integrated <- FindNeighbors(seu.integrated,dims = c(1:20))
seu.integrated <- FindClusters(seu.integrated,resolution = 0.2)

####-------5.2 assign clusters-------------

TE.markers <- c("Cdx2","Krt8","Krt18","Ascl2","Tacstd2")
ICM.markers <- c("Pou5f1","Nanog","Sox2","Esrrb","Sox15")
PE.markers <- c("Gata4","Gata6","Sox17","Pdgfra","Col4a1")
p1 <- UMAPPlot(seu.integrated) 
p1
tmp.markers <- c(TE.markers,
                 "Gapdh",ICM.markers,
                 "Ppia",PE.markers)

# p2 <- FeaturePlot(seu.integrated,features = c(TE.markers,
#                                    "Gapdh",ICM.markers,
#                                    "Ppia",PE.markers),combine = F) 

p2 <- lapply(tmp.markers, function(x){
  res.p <- FeaturePlot(seu.integrated,features = x)+
    scale_color_gradientn(colours = rdwhbu(256))
  return(res.p)
})

p <- wrap_plots(c(list(p1),p2),nrow = 3) & 
  theme_cowplot(font_size = 24) &
  NoAxes()

ggsave(plot = p,filename = myFileName(prefix = "res/fig/Fig3_EPS_blastoid_markers_raw",
                             suffix = ".png"),
       dpi = 350,width = 18,height = 8)

UMAPPlot(seu.integrated,label=T)
?UMAPPlot
seu.integrated <- RenameIdents(seu.integrated,
                               "0"="intermediate",
                               "1"="PE",
                               "2"="intermediate",
                               "3"="PE",
                               "4"="intermediate",
                               "5"="EPI",
                               "6"="TE",
                               "7"="EPI")
levels(seu.integrated) <- c("EPI","TE","PE","intermediate")
p1 <- UMAPPlot(seu.integrated) + 
  scale_color_manual(values = c("#e90b8e","#3b53a5","#6ebe45","grey"))
p1


p <- wrap_plots(c(list(p1),p2),nrow = 3) & 
  theme_cowplot(font_size = 24) &
  NoAxes()
p
ggsave(plot = p,filename = myFileName(prefix = "res/fig/Fig3_EPS_blastoid_markers",
                                      suffix = ".png"),
       dpi = 350,width = 18,height = 8)


p3 <- UMAPPlot(seu.integrated,pt.size=1.5,
               split.by="orig.ident") +
  scale_color_manual(values = c("#e90b8e","#3b53a5","#6ebe45","grey"))+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
p3
ggsave(filename = myFileName(prefix = "res/fig/Fig3_EPS_split_plot",
                             suffix = ".png"),
       width = 8,height = 4,dpi = 350)

saveRDS(seu.integrated,file = myFileName(prefix = "res/R/EPS_blastoid",suffix = ".rds"))


#####----show markers-----
TE.markers <- c("Cdx2","Krt8","Krt18","Ascl2","Tacstd2")
ICM.markers <- c("Pou5f1","Nanog","Sox2","Esrrb","Sox15")
PE.markers <- c("Gata4","Gata6","Sox17","Pdgfra","Col4a1")
tmp.markers <- c("Gapdh","Ppia",
                 TE.markers,
                 ICM.markers,
                 PE.markers)
tmp.seu <- seu.integrated
DefaultAssay(tmp.seu) <- "RNA"
####myscAreaPlot(tmp.seu = tmp.seu,features = "Sox2",tmp.color = divergentcolor(length(levels(tmp.seu))))
levels(tmp.seu) <- c("TE","EPI","PE")
myscAreaPlot(tmp.seu = tmp.seu,
             features = tmp.markers,
             tmp.color = divergentcolor(length(levels(tmp.seu))))

###------------5.3 IS -----------

###---------5.3.1 define function-------
MyAverageSeurat <- function(seu.obj){
  cat("Calcualte Average Expression",sep = "\n")
  tmp.exp.data  <-  GetAssayData(seu.obj,slot = "data",assay = "RNA")
  tmp.exp.meta  <- seu.obj[[]]
  tmp.cell.type <- levels(seu.obj)
  tmp.res <- lapply(tmp.cell.type,FUN = function(x){
    cat(x,sep = "\n")
    tmp.id <- tmp.exp.meta %>% 
      dplyr::filter(cell_type == x) %>%
      pull(cell_id)
    tmp.exp.data <- rowMeans(expm1(tmp.exp.data[,tmp.id]))
    return(tmp.exp.data)
  })
  tmp.res <- Reduce("cbind",tmp.res)
  colnames(tmp.res) <- tmp.cell.type
  return(tmp.res)
}


myInteraction_Score <- function(gene.exp.mat=tmp.mat,
                                cell.type=tmp.cell.type,
                                Lgenelist=Lgenelist,
                                Rgenelist=Rgenelist){
  cat("Calcualte IS",sep = "\n")
  tmp.res.IS <- lapply(cell.type,FUN = function(x){
    cat(x,sep = "\n")
    res <- gene.exp.mat[Lgenelist,x] * gene.exp.mat[Rgenelist,]
    colnames(res) <- paste0(x,"_",cell.type)
    return(res)
  })
  tmp.res.IS <- Reduce(cbind,tmp.res.IS)
  rownames(tmp.res.IS) <- paste0(Lgenelist,"_",Rgenelist)
  return(tmp.res.IS)
}

myRandomIS <- function(seu.obj=seu,Lgenelist=Lgenelist,Rgenelist=Rgenelist){
  cat("random_sample",sep = "\n")
  tmp.exp.data  <-  GetAssayData(seu.obj,slot = "data",assay = "RNA")
  tmp.cell.type <- levels(seu.obj)
  tmp.label <- unname(seu.obj$cell_type)
  tmp.label <- sample(x = as.character(tmp.label),size = length(tmp.label),replace = F)
  seu.obj$random_label <- tmp.label
  tmp.exp.meta  <- seu.obj[[]]
  
  tmp.mat <- lapply(tmp.cell.type,FUN = function(x){
    cat(x,sep = "\n")
    tmp.id <- tmp.exp.meta %>% 
      dplyr::filter(random_label == x) %>%
      pull(cell_id)
    tmp.exp.data <- rowMeans(expm1(tmp.exp.data[,tmp.id]))
    return(tmp.exp.data)
  })
  tmp.mat <- Reduce("cbind",tmp.mat)
  colnames(tmp.mat) <- tmp.cell.type
  
  tmp.res <- myInteraction_Score(gene.exp.mat = tmp.mat,
                                 cell.type = colnames(tmp.mat),
                                 Lgenelist = Lgenelist,
                                 Rgenelist = Rgenelist)
  return(tmp.res)
}

myISwithPvalue <- function(seu=seu,Lgenelist = Lgenelist,Rgenelist = Rgenelist,calc.p=FALSE){
  
  tmp.seu <- seu
  tmp.res <- MyAverageSeurat(seu.obj = tmp.seu)
  tmp.res.IS <- myInteraction_Score(gene.exp.mat = tmp.res,
                                    cell.type = colnames(tmp.res),
                                    Lgenelist = Lgenelist,
                                    Rgenelist = Rgenelist)
  tmp.res.list <- list(IS=tmp.res.IS)
  
  if(calc.p){
    #### use future to do parallel
    options(future.globals.maxSize= 1024^3)
    plan(multisession, workers = 10)
    nboot = 1000
    tmp.res.IS.boot <- future_replicate(nboot, myRandomIS(seu.obj = seu,Lgenelist = Lgenelist,Rgenelist = Rgenelist))
    .tempGetPvalueFromOnePairs <- function(xx,tmp.res.IS,tmp.res.IS.boot,nboot=1000){
      test.p <- lapply(1:nrow(tmp.res.IS), 
                       FUN = function(x) sum(tmp.res.IS[x,xx] >= tmp.res.IS.boot[x,xx,])/nboot)
      test.p <- unlist(test.p)
      return(test.p)
    }
    GetPvalueFromISBoot <- function(tmp.res.IS,tmp.res.IS.boot){
      tmp.res <- lapply(1:ncol(tmp.res.IS),
                        FUN = function(x) .tempGetPvalueFromOnePairs(xx = x,tmp.res.IS = tmp.res.IS,tmp.res.IS.boot = tmp.res.IS.boot))
      tmp.res.df <- Reduce(cbind,tmp.res)
      return(tmp.res.df)
    }
    tmp.res.IS.Pvalue <- GetPvalueFromISBoot(tmp.res.IS = tmp.res.IS,tmp.res.IS.boot = tmp.res.IS.boot)
    rownames(tmp.res.IS.Pvalue) <- rownames(tmp.res.IS)
    colnames(tmp.res.IS.Pvalue) <- colnames(tmp.res.IS)
  }
  tmp.res.list <- c(tmp.res.list,list(IS.p=tmp.res.IS.Pvalue))
  return(tmp.res.list)
}


####--------5.3.2 load data----------

####load seurat
seu.integrated <- readRDS(file = "res/R/EPS_blastoid_20211017.rds")
seu.integrated$cell_type = Idents(seu.integrated)
seu.integrated$cell_id = colnames(seu.integrated)
seu.integrated <- subset(seu.integrated,idents = "intermediate", invert = TRUE)

seu.list <- SplitObject(seu.integrated,split.by = 'orig.ident')

###load pairs
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)



####filter genes
DefaultAssay(seu.integrated) <- "RNA"
gene_symbols <- rownames(seu.integrated)
l.remove <- setdiff(Lgenelist,gene_symbols)
r.remove <- setdiff(Rgenelist,gene_symbols)
index.remove <- c(which(Lgenelist %in% l.remove),which(Rgenelist %in% r.remove))
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

####load cell.type
# tmp.cell.type <- c("EPI","TE","PE")
# 
# 
# head(colnames(seu.integrated))
# tmp.seu <- seu.list$EPS_blastoid
# tmp.res <- MyAverageSeurat(seu.obj = tmp.seu)
# tmp.res.IS <- myInteraction_Score(gene.exp.mat = tmp.res[,tmp.cell.type],
#                                   cell.type = tmp.cell.type,
#                                   Lgenelist = Lgenelist,
#                                   Rgenelist = Rgenelist)
# rownames(tmp.res.IS) <- paste0(Lgenelist,"_",Rgenelist)


tmp.seu <- seu.list$EPS_blastoid
cell_id <- sample(x = colnames(tmp.seu),size = ncol(seu.list$blastocyst))
tmp.seu <- subset(tmp.seu,cells=cell_id)
UMAPPlot(tmp.seu)
UMAPPlot(seu.list$blastocyst)
tmp.seu 

tmp.res.list <- myISwithPvalue(seu = tmp.seu,
                               Lgenelist = Lgenelist,
                               Rgenelist = Rgenelist,
                               calc.p=T)

tmp.res.list.EPS_blastoid <- tmp.res.list

tmp.res.list <- myISwithPvalue(seu = seu.list$blastocyst,
                               Lgenelist = Lgenelist,
                               Rgenelist = Rgenelist,
                               calc.p=T)
tmp.res.list.blastocyst <- tmp.res.list


####------5.3.3 buble plot--------
head(tmp.res.list$IS)
head(tmp.res.list$IS.p)

tmp.IS.p.buble.plot <- function(tmp.res.list,title=NULL){
  
  data.plot <- as.data.frame(tmp.res.list$IS)
  data.plot <- data.plot %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS",-LRpairs)
  
  data.plot.p <- as.data.frame(tmp.res.list$IS.p)
  data.plot.p <- data.plot.p %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS.p",-LRpairs)
  
  
  data.plot$IS.p <- data.plot.p$IS.p
  
  #### load eLR
  
  tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
  tmp.cutoff <- 0.1
  tmp.LR.pairs.1 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC > tmp.cutoff) %>%
    head(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs.2 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC < tmp.cutoff) %>%
    tail(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs <- c(tmp.LR.pairs.1,tmp.LR.pairs.2)
  
  tmp.LR.pairs <- gsub(pattern = "_",replacement = "-",tmp.LR.pairs)
  data.plot$LRpairs <- gsub(pattern = "_",replacement = "-",data.plot$LRpairs)
  data.plot$cell_pairs <- gsub(pattern = "_",replacement = "-",data.plot$cell_pairs)
  
  #tmp.LR.pairs <- c(tmp.LR.pairs,"Bmp4_Bmpr2")
  
  # tmp.LR.pairs.1
  # tmp.LR.pairs.2
  
  plot.data <- data.plot %>%
    filter(LRpairs %in% tmp.LR.pairs)
  plot.data$LRpairs <- factor(plot.data$LRpairs,levels = rev(tmp.LR.pairs))
  
  ###my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  tmp.min <- 0.001
  #log10(tmp.min)
  plot.data$IS.p[plot.data$IS.p <= tmp.min] = tmp.min
  
  tmp.min <- 0.000001
  plot.data$IS[plot.data$IS <= tmp.min] = tmp.min
  
  
  p <- ggplot(plot.data,aes(x=cell_pairs,y=LRpairs)) +
    geom_point(aes(size=-log10(IS.p),color=log10(IS))) +
    scale_color_gradientn('Log10(IS)', colors=rdbu(100)) +
    theme_bw(base_size = 24) +
    ggtitle(title)+
    guides(size = guide_legend(order = 1),
           color = guide_colorbar(order = 2))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=24, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title=element_blank(),plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  return(p)
  
}




p1 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.blastocyst,title = "blastocyst")
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_blastocyst_buble",suffix = ".png"),
       width = 8,height = 16,dpi = 350)

p2 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.EPS_blastoid,title = "EPS blastoid")
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_EP_buble",suffix = ".png"),
       width = 8,height = 16,dpi = 350)

p1+theme(legend.position = "none",
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18,face = "italic")) | 
  p2+theme(axis.text.x = element_text(size = 18),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank())

ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS",suffix = ".png"),
       width = 9,height = 12,dpi = 350,units = "in")


####------5.3.4 diff---------
tmp.res.list <- tmp.res.list.blastocyst
data.plot <- as.data.frame(tmp.res.list$IS)
data.plot <- data.plot %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "cell_pairs",value="IS",-LRpairs)

data.plot.p <- as.data.frame(tmp.res.list$IS.p)
data.plot.p <- data.plot.p %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "cell_pairs",value="IS.p",-LRpairs)


data.plot$IS.p <- data.plot.p$IS.p

#### load eLR

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
tmp.cutoff <- 0.1
tmp.LR.pairs.1 <- tmp.res.cor.1 %>%
  dplyr::filter(PCC > tmp.cutoff) %>%
  pull(LRpairs)

tmp.LR.pairs.2 <- tmp.res.cor.1 %>%
  dplyr::filter(PCC < tmp.cutoff) %>%
  pull(LRpairs)

tmp.LR.pairs <- c(tmp.LR.pairs.1,tmp.LR.pairs.2)

#tmp.LR.pairs <- c(tmp.LR.pairs,"Bmp4_Bmpr2")

tmp.LR.pairs.1
tmp.LR.pairs.2

plot.data <- data.plot %>%
  filter(LRpairs %in% tmp.LR.pairs)



tmp.plot <- log10(tmp.res.list.EPS_blastoid$IS+1) - log10(tmp.res.list.blastocyst$IS+1)
tmp.LR.pairs <- intersect(rownames(tmp.plot),tmp.LR.pairs)
tmp.plot <- tmp.plot[tmp.LR.pairs,]
tmp.sd <- -sort(-apply(tmp.plot, 1, sd))
tmp.sd <- tmp.sd[tmp.sd>0.1]

tmp.plot <- tmp.plot[names(tmp.sd),]
tmp.plot[tmp.plot>1] <- 1
tmp.plot[tmp.plot < -1] <- -1

png(filename = myFileName(prefix = "res/fig/fig3_scInt_diff",suffix = ".png"),
    width = 8,height = 16,units = "in",res = 350)
ComplexHeatmap::pheatmap(tmp.plot,
         color = rdbu(100),
         cluster_cols = F,
         show_rownames = T,fontsize = 20,
         heatmap_legend_param = list(title = expression(Delta(IS)),
                                     title_gp = gpar(col = "black", fontsize = 20),
                                     labels_gp = gpar(col = "black", fontsize = 20)),
         show_colnames = T)
dev.off()




#####----------5.3.5 diff pheatmap---------------
tmp.IS.bubble.data.get <- function(tmp.res.list,title=NULL){
  
  data.plot <- as.data.frame(tmp.res.list$IS)
  data.plot <- data.plot %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS",-LRpairs)
  
  data.plot.p <- as.data.frame(tmp.res.list$IS.p)
  data.plot.p <- data.plot.p %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS.p",-LRpairs)
  
  
  data.plot$IS.p <- data.plot.p$IS.p
  
  #### load eLR
  
  tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
  tmp.cutoff <- 0.1
  tmp.LR.pairs.1 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC > tmp.cutoff) %>%
    head(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs.2 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC < tmp.cutoff) %>%
    tail(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs <- c(tmp.LR.pairs.1,tmp.LR.pairs.2)
  
  tmp.LR.pairs <- gsub(pattern = "_",replacement = "-",tmp.LR.pairs)
  data.plot$LRpairs <- gsub(pattern = "_",replacement = "-",data.plot$LRpairs)
  data.plot$cell_pairs <- gsub(pattern = "_",replacement = "-",data.plot$cell_pairs)
  
  #tmp.LR.pairs <- c(tmp.LR.pairs,"Bmp4_Bmpr2")
  
  # tmp.LR.pairs.1
  # tmp.LR.pairs.2
  
  plot.data <- data.plot %>%
    filter(LRpairs %in% tmp.LR.pairs)
  plot.data$LRpairs <- factor(plot.data$LRpairs,levels = rev(tmp.LR.pairs))
  
  ###my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  tmp.min <- 0.001
  #log10(tmp.min)
  plot.data$IS.p[plot.data$IS.p <= tmp.min] = tmp.min
  
  tmp.min <- 0.000001
  plot.data$IS[plot.data$IS <= tmp.min] = tmp.min
  
  return(plot.data)
  
}

p1 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.blastocyst,title = "blastocyst")
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_blastocyst_buble",suffix = ".png"),
       width = 8,height = 16,dpi = 350)

p2 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.EPS_blastoid,title = "EPS blastoid")
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_EP_buble",suffix = ".png"),
       width = 8,height = 16,dpi = 350)

tmp.data.get.1 <- tmp.IS.bubble.data.get(tmp.res.list = tmp.res.list.blastocyst)
tmp.data.get.2 <- tmp.IS.bubble.data.get(tmp.res.list = tmp.res.list.EPS_blastoid)
all(tmp.data.get.1$LRpairs==tmp.data.get.2$LRpairs)
tmp.plot <- tmp.data.get.1 %>%
  mutate(IS.diff=tmp.data.get.1$IS-tmp.data.get.2$IS) %>%
  dplyr::select(LRpairs,cell_pairs,IS.diff) %>%
  spread(key = cell_pairs,value = IS.diff) %>%
  column_to_rownames("LRpairs") 

tmp.plot <- tmp.plot[rev(1:nrow(tmp.plot)),]
tmp.plot[tmp.plot>1] <- 1
tmp.plot[tmp.plot < -1] <- -1

png(filename = myFileName(prefix = "res/fig/fig3_scInt_diff",suffix = ".png"),
    width = 8,height = 16,units = "in",res = 350)
ComplexHeatmap::Heatmap(matrix = tmp.plot,
                        col = rdbu(100),
                        cluster_columns = F,
                        cluster_rows = F,
                        row_names_gp = gpar(fontsize=20,fontface="italic"),
                        column_names_gp = gpar(fontsize=20),
                        rect_gp = gpar(col = "black"),
                        heatmap_legend_param = list(title = expression(Delta(IS)),
                                                    title_gp = gpar(col = "black", fontsize = 20),
                                                    labels_gp = gpar(col = "black",fontsize = 20)),
                        show_column_names = T,
                        show_row_names = T)
dev.off()
?Heatmap
?expression

####----5.3.6 circos plot--------


tmp.res.list <- tmp.res.list.blastocyst

tmp.circos.plot <- function(tmp.res.list){
  data.plot <- as.data.frame(tmp.res.list$IS)
  data.plot <- data.plot %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS",-LRpairs)
  
  data.plot.p <- as.data.frame(tmp.res.list$IS.p)
  data.plot.p <- data.plot.p %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS.p",-LRpairs)
  data.plot$IS.p <- data.plot.p$IS.p
  
  data.plot$Ligand <- unlist(lapply(strsplit(data.plot$LRpairs,split = "_"),FUN = function(x) x[1])) 
  data.plot$Receptor <- unlist(lapply(strsplit(data.plot$LRpairs,split = "_"),FUN = function(x) x[2])) 
  data.plot$cell_L <- unlist(lapply(strsplit(data.plot$cell_pairs,split = "_"),FUN = function(x) x[1])) 
  data.plot$cell_R <- unlist(lapply(strsplit(data.plot$cell_pairs,split = "_"),FUN = function(x) x[2]))
  
  
  data.plot.gene <- data.frame(cell=c(data.plot$cell_L,data.plot.circos$cell_R),
                               gene=c(data.plot$Ligand,data.plot.circos$Receptor),
                               lr=c(rep("ligand",length(data.plot$Ligand)),
                                    rep("receptor",length(data.plot$Receptor))),
                               stringsAsFactors = F) %>%
    unique()
  data.plot.gene$gene_id <- paste0("gene",1:nrow(data.plot.gene))
  data.plot.gene$map_id <- paste0(data.plot.gene$cell,"_",data.plot.gene$gene,"_",data.plot.gene$lr)
  
  
  
  library(circlize)
  df <- data.plot.gene %>%
    arrange(cell)
  
  tmp.cell.type <- unique(df$cell)
  tmp.cell.colors <- c("#e90b8e","#3b53a5","#6ebe45","grey")
  names(tmp.cell.colors) <- c("EPI","TE","PE","intermediate")
  
  circos.clear()
  circos.par(canvas.xlim =c(-1.1,1.1),
             canvas.ylim = c(-1.1,1.1),
             gap.degree=0.001,
             cell.padding = c(0.02, 0, 0.02, 0))
  fa = df$gene_id
  fa = factor(fa,levels = fa)
  circos.initialize(factors = fa, xlim = c(0,1))
  
  circos.trackPlotRegion(
    
    ylim = c(0, 1), track.height = 0.15, bg.border = NA,
    
    panel.fun = function(x, y) {
      
      sector.index = get.cell.meta.data("sector.index")
      
      xlim = get.cell.meta.data("xlim")
      
      ylim = get.cell.meta.data("ylim")
      
    } )
  
  
  
  lapply(1:length(tmp.cell.type),FUN = function(ii){
    x <- tmp.cell.type[ii]
    tmp.gene.id <- df %>%
      dplyr::filter(cell==x) %>%
      pull(gene_id)
    highlight.sector(tmp.gene.id, 
                     track.index = 1,
                     text = x, 
                     niceFacing = F, 
                     font = 2, 
                     col = tmp.cell.colors[x])
  })
  
  circos.trackPlotRegion(
    
    ylim = c(0, 1), track.height = 0.08, bg.border = NA,
    
    panel.fun = function(x, y) {
      
      sector.index = get.cell.meta.data("sector.index")
      
      xlim = get.cell.meta.data("xlim")
      
      ylim = get.cell.meta.data("ylim")
      
    } )
  
  
  
  tmp.cat <- unique(df$lr)
  tmp.cat.color <- c('#ffc000',"#00af50")
  names(tmp.cat.color) <- c("ligand","receptor")
  lapply(1:length(tmp.cat),FUN = function(ii){
    x <- tmp.cat[ii]
    tmp.gene.id <- df %>%
      dplyr::filter(lr==x) %>%
      pull(gene_id)
    highlight.sector(tmp.gene.id, 
                     track.index = 2,
                     text = x, niceFacing = F,
                     col = tmp.cat.color[x],
                     text.col = 'black')
  })
  
  
  lrid <- data.plot %>%
    dplyr::filter(IS.p < 0.01) %>%
    dplyr::mutate(xx = paste0(cell_L,"_",Ligand,"_","ligand"),
                  yy = paste0(cell_R,"_",Receptor,"_","receptor")) %>%
    dplyr::mutate(xx=plyr::mapvalues(xx,
                                     from = data.plot.gene$map_id,
                                     to = data.plot.gene$gene_id),
                  yy=plyr::mapvalues(yy,
                                     from = data.plot.gene$map_id,
                                     to = data.plot.gene$gene_id))
  
  
  
  for(i in 1:nrow(lrid)){
    circos.link(sector.index1 = lrid$xx[i], 
                point1 = 0, 
                sector.index2 = lrid$yy[i], 
                point2 = 0 ,
                directional = 0,
                h=0.5, 
                lwd=0.1, 
                col="grey",
                lty=1)
  }
  
}

tmp.circos.plot(tmp.res.list = tmp.res.list.blastocyst)

p <- recordPlot()

png(filename = myFileName(prefix = "res/fig/fig3_circos_blastocyst",suffix = ".png"),width = 8,height = 8,units = "in",res = 350)
p
title("blastocyst")
dev.off()


tmp.circos.plot(tmp.res.list = tmp.res.list.EPS_blastoid)
p <- recordPlot()
png(filename = myFileName(prefix = "res/fig/fig3_circos_EPS_blastoid",suffix = ".png"),width = 8,height = 8,units = "in",res = 350)
p
title("EPS_blastoid")
dev.off()



####-------6. validate on ZGA gene---------------
####20210601

#### load expression data

tmp.df.1 <- read.delim(file = "data/GSE71434_FPKM_inhibition.txt",sep = "\t",stringsAsFactors = F)
tmp.df.2 <- read.delim(file = "data/GSE71434_FPKM_stage.txt",sep = "\t",stringsAsFactors = F)
tmp.df <- merge(tmp.df.1,tmp.df.2,by="Gene")


#### load eLR data
tmp.cutoff <- 0.1
LRpairs.df <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
LRpairs.df <- LRpairs.df %>%
  dplyr::filter(abs(PCC) >= tmp.cutoff)

L.gene.list <- LRpairs.df$Lgene
R.gene.list <- LRpairs.df$Rgene

tmp.df <- tmp.df[!duplicated(tmp.df$Gene),]
data.plot <-  tmp.df %>% 
  dplyr::select(Gene,
                X2cell_early,
                X2cell_late,
                X2cell_alpha_amanitin) 
rownames(data.plot) <- data.plot$Gene
data.plot <- data.plot[,-1]
colnames(data.plot) <- c("early_2cell","late_2cell","2cell_alpha_amanitin")
#data.plot <- log10(data.plot+1)

data.plot.L <- data.plot[LRpairs.df$Lgene,]
data.plot.R <- data.plot[LRpairs.df$Rgene,]

data.plot.IS <- data.plot.L * data.plot.R
rownames(data.plot.IS) <- paste0(LRpairs.df$Lgene,"_",LRpairs.df$Rgene)

colnames(data.plot.IS)
data.plot.IS <- log10(data.plot.IS+1)

num_clusters = 6
ph <- pheatmap::pheatmap(data.plot.IS,
               scale = "none",
               color = rdbu(100),
               show_rownames = F,
               cutree_rows = num_clusters,
               silent = T)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,num_clusters)))

png(filename = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_pheatmap",suffix = ".png"),
    width = 6,height = 12,units = "in",res = 350)
pheatmap::pheatmap(data.plot.IS,
         scale = "none",
         color = rdbu(100),
         show_rownames = F,
         cutree_rows = num_clusters,
         angle_col = "45",
         fontsize = 28,
         annotation_row = annotation_row)
dev.off()


png(filename = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_pheatmap_filter",suffix = ".png"),
    width = 6,height = 8,units = "in",res = 350)
tmp.idx <- annotation_row %>%
  dplyr::filter(Cluster %in% c(4,6)) %>%
  rownames()
num_clusters <- 4
tmp.data.plot.ttt <- data.plot.IS[tmp.idx,]
#tmp.data.plot.ttt[tmp.data.plot.ttt>2] <- 2
ph <- pheatmap::pheatmap(tmp.data.plot.ttt,
               scale = "none",
               color = rdbu(100),
               show_rownames = T,
               show_colnames = T,
               cutree_rows = num_clusters,
               annotation_row = annotation_row,
               silent = T,
               angle_col = 45)
annotation_row.1 <- data.frame(Cluster = factor(cutree(ph$tree_row,num_clusters)))
ph <- pheatmap::pheatmap(tmp.data.plot.ttt,
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = F,
                         show_colnames = T,
                         cutree_rows = num_clusters,
                         annotation_row = annotation_row.1,
                         silent = T,
                         fontsize = 28,
                         angle_col = 45)
ph
dev.off()


tmp.idx <- annotation_row.1 %>%
  dplyr::filter(Cluster %in% c(2,3,4)) %>%
  rownames()

png(filename = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".png"),width = 3,height = 9,units = "in",res = 350)
pheatmap::pheatmap(data.plot.IS[tmp.idx,],
         scale = "none",
         color = rdbu(100),
         show_rownames = T,
         show_colnames = T,
         silent = F,
         angle_col = 45)
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".pdf"),
    width = 6,height = 12)
ComplexHeatmap::pheatmap(data.plot.IS[tmp.idx,],
                   scale = "none",
                   color = rdbu(100),
                   show_rownames = T,
                   show_colnames = T,
                   fontsize = 16,
                   name = "IS",
                   angle_col = "45")
dev.off()


ph <- pheatmap::pheatmap(data.plot.IS[tmp.idx,],
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = T,
                         show_colnames = T,
                         cellwidth = 16,
                         cellheight = 16,
                         fontsize = 16,
                         angle_col = "45")

get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}
tmp.size <- get_plot_dims(ph)

png(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".png"),
    width = tmp.size$width+0.5,
    height = tmp.size$height+0.5,
    units = "in",res = 350)
ComplexHeatmap::pheatmap(data.plot.IS[tmp.idx,],
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = T,
                         show_colnames = T,
                         fontsize = 16,
                         heatmap_legend_param = list(title_gp=gpar(fontsize=16,fontface="bold"),
                                                     labels_gp=gpar(fontsize=16,fontface="bold")),
                         name = "IS",
                         angle_col = "45")
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".pdf"),
    width = tmp.size$width+0.5,
    height = tmp.size$height+0.5)
ComplexHeatmap::pheatmap(data.plot.IS[tmp.idx,],
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = T,
                         show_colnames = T,
                         fontsize = 16,
                         heatmap_legend_param = list(title_gp=gpar(fontsize=16,fontface="bold"),
                                                     labels_gp=gpar(fontsize=16,fontface="bold")),
                         name = "IS",
                         angle_col = "45")
dev.off()


tmp.data.plot <- data.plot.IS[tmp.idx,]
rownames(tmp.data.plot) <- gsub(rownames(tmp.data.plot),pattern = "_",replacement = "-")
colnames(tmp.data.plot) <- c("early 2-cell","late 2-cell","2-cell alpha amanitin")

tmp.data.plot[tmp.data.plot > 2] <- 2
pdf(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".pdf"),
    width = tmp.size$width+2.5,
    height = tmp.size$height+0.5)
ComplexHeatmap::Heatmap(matrix = tmp.data.plot,
                        col = rdbu(100),
                        row_names_gp = gpar(fontsize=20,col=divergentcolor(nrow(tmp.data.plot)),
                                            fontface="italic"),
                        rect_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize=20),
                        column_names_rot = 45,
                        column_names_centered = F,
                        heatmap_legend_param = list(title="IS",
                                                    title_gp=gpar(fontsize=20,fontface="bold"),
                                                    labels_gp=gpar(fontsize=16,fontface="bold")))
dev.off()


png(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".png"),
    width = tmp.size$width+2,
    height = tmp.size$height+2,
    units = "in",res = 350)
ComplexHeatmap::Heatmap(matrix = tmp.data.plot,
                        col = rdbu(100),
                        row_names_gp = gpar(fontsize=20,
                                            col=ifelse(rownames(tmp.data.plot) %in% c("Fgf4-Fgfr1","Fgf4-Fgfr2"),"red","black"),
                                            fontface="italic"),
                        rect_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize=20),
                        column_names_rot = 45,
                        column_names_centered = F,
                        heatmap_legend_param = list(title="IS",
                                                    title_gp=gpar(fontsize=20,fontface="bold"),
                                                    labels_gp=gpar(fontsize=20)))
dev.off()
