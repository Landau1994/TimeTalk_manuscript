
myCellChatToSeu <- function(object,slot.use="data"){
  meta <- object@meta
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  }
  else {
    meta$group.cellchat <- object@idents
  }
  w10x <- Seurat::CreateSeuratObject(counts = slot(object = object,name = slot.use),
                                     meta.data = meta)
  Seurat::Idents(w10x) <- "group.cellchat"
  return(w10x)
}

myscAreaPlot <- function(tmp.seu = tmp.seu, 
                         features = "Col1a1",
                         tmp.color= tmp.color){
  data <- Seurat::FetchData(object = tmp.seu, vars = features, slot = "data")
  data$label <- Idents(tmp.seu)
  data <- data %>%
    gather(key = "features",value = "gene.exp",-label)
  
  data <- data %>%
    group_by(features,label) %>%
    mutate(feature = (gene.exp-min(gene.exp))/(max(gene.exp)-min(gene.exp))) %>%
    arrange(-feature) %>%
    mutate(index = row_number()) %>%
    mutate(index = (index - min(index))/(max(index)-min(index))) %>%
    ungroup()
  
  data$features <- factor(data$features,levels = features)
  
  p <- ggplot(data = data,aes(x=-index,y= feature ,fill=label))+
    geom_area(color=NA)+
    facet_grid(features~label)+
    scale_fill_manual(values = tmp.color)+
    xlab(NULL)+
    ylab("Relative gene expression level per cluster")+
    scale_y_continuous(expand = c(0,0))+
    theme_few()+
    theme(axis.line  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.background.x = element_rect(color = "black",fill="white"),
          strip.background.y = element_rect(color = "black",fill="white"),
          strip.text.y = element_text(angle = 0,face = "italic"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing = unit(0,"lines"),
          legend.position = "none")
  p
  return(p)
}



myPlotCellChatFeature <- function (object, features = NULL, signaling = NULL, enriched.only = TRUE, 
                                   using.slot="data",
                                   type = c("violin", "dot","scArea"), color.use = NULL, 
                                   group.by = NULL, ...) 
{
  type <- match.arg(type)
  meta <- object@meta
  
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  }
  else {
    meta$group.cellchat <- object@idents
  }
  w10x <- Seurat::CreateSeuratObject(counts = slot(object = object,name = using.slot), 
                                     meta.data = meta)
  if (is.null(group.by)) {
    group.by <- "group.cellchat"
  }
  Seurat::Idents(w10x) <- group.by
  
  if (is.null(color.use)) {
    numCluster <- length(levels(Seurat::Idents(w10x)))
    color.use <- scPalette(numCluster)
  }
  if (!is.null(features) & !is.null(signaling)) {
    warning("`features` will be used when inputing both `features` and `signaling`!")
  }
  if (!is.null(features)) {
    feature.use <- features
  }
  else if (!is.null(signaling)) {
    res <- extractEnrichedLR(object, signaling = signaling, 
                             geneLR.return = TRUE, enriched.only = enriched.only)
    feature.use <- res$geneLR
  }
  if (type == "violin") {
    gg <- StackedVlnPlot(w10x, features = feature.use, color.use = color.use, 
                         ...)
  }
  if (type == "dot") {
    gg <- dotPlot(w10x, features = feature.use, ...)
  }
  if (type == "scArea"){
    gg <- myscAreaPlot(w10x,features = feature.use,tmp.color = color.use)
  }
  return(gg)
}



####show color bar

my_showcolorbar <- function(tmp.colorbar,tmp.name,n=100){
  image(
    1:n, 1, as.matrix(1:n),
    col = tmp.colorbar(n),
    xlab = tmp.name, 
    ylab = "", 
    xaxt = "n", 
    yaxt = "n", 
    bty = "n"
  )
}


####Define function


my_getLR <- function(tmp.L.exp = data.plot[Lgene.list,tmp.id],
         tmp.R.exp = data.plot[Rgene.list,tmp.id],
         tmp.p.cutoff=0.05,tmp.cell.cutoff=3){
  
  cell_number <- ncol(tmp.L.exp)
  
  tmp.res.PCC <- sapply(1:nrow(tmp.L.exp), function(i) {
    cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "pearson")
  })
  tmp.res.PCC.p <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp <- cor.test(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "pearson")
    return(tmp$p.value)
  })
  
  tmp.res.SCC <- sapply(1:nrow(tmp.L.exp), function(i) {
    cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "spearman")
  })
  tmp.res.SCC.p <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp <- cor.test(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "spearman")
    return(tmp$p.value)
  })
  
  tmp.res.KCC <- sapply(1:nrow(tmp.L.exp), function(i) {
    cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "kendall")
  })
  tmp.res.KCC.p <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp <- cor.test(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "kendall")
    return(tmp$p.value)
  })
  
  tmp.res.rate.L <- sapply(1:nrow(tmp.L.exp), function(i) sum(as.numeric(tmp.L.exp[i,])>0)/length(as.numeric(tmp.L.exp[i,])>0))
  tmp.res.rate.R <- sapply(1:nrow(tmp.R.exp), function(i) sum(as.numeric(tmp.R.exp[i,])>0)/length(as.numeric(tmp.R.exp[i,])>0))
  tmp.res.cor.1 <- data.frame(Lgene=Lgene.list,
                              Rgene=Rgene.list,
                              LRpairs=paste0(Lgene.list,"_",Rgene.list),
                              PCC=tmp.res.PCC,
                              PCC.p = tmp.res.PCC.p,
                              SCC = tmp.res.SCC,
                              SCC.p = tmp.res.SCC.p,
                              KCC = tmp.res.KCC,
                              KCC.p = tmp.res.KCC.p,
                              detection.rate.L=tmp.res.rate.L,
                              detection.rate.R=tmp.res.rate.R,
                              stringsAsFactors = F) %>%
    dplyr::filter(detection.rate.L > tmp.cell.cutoff/cell_number & detection.rate.R > tmp.cell.cutoff/cell_number)
  
  tmp.res.cor.1$PCC[is.na(tmp.res.cor.1$PCC)] <- 0 
  tmp.res.cor.1 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC.p < tmp.p.cutoff) %>%
    arrange(-PCC) 
  
  return(tmp.res.cor.1)
  
}




####define function
myenrichGO_plot <- function(tmp.res.plot,
                            p.adjust.cut=0.01,
                            show_number=20,
                            title="The Most Enriched GO Terms",
                            fill.color="#8DA1CB",
                            plot.ylab="BP",
                            font.size = 12){
  tmp.res.plot$logp.adjust <- -log10(tmp.res.plot$p.adjust)
  tmp.res.plot <- tmp.res.plot%>%
    dplyr::filter(p.adjust < p.adjust.cut) %>%
    slice_max(order_by = logp.adjust,
              n=show_number,
              with_ties=F)
  tmp.res.plot$Description <- factor(tmp.res.plot$Description,
                                     levels = rev(tmp.res.plot$Description))
  p <- ggplot() + 
    geom_bar(data = tmp.res.plot, aes_string(x = "Description", 
                                             y = "logp.adjust"),
             stat = "identity",width=0.5,
             fill=fill.color)+
    coord_flip() +
    ggtitle(title)+
    xlab(plot.ylab)+
    ylab("-log10(p.adjust)")+
    theme_cowplot(font_size = font.size) +
    theme(axis.ticks.y = element_blank(),
          plot.title = element_text(hjust  = 0.5))
  
  return(p)
}



my_LRScatterPlot <- function(mat = data.plot,
                          mat.meta = tmp.meta,
                          tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                          tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                          PCC = tmp.res.cor.1$PCC[tmp.idx]){
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(Lgene = mat.L,
                       Rgene = mat.R,
                       stringsAsFactors = F)
  colnames(res.df) <- c(tmp.Lgene,tmp.Rgene)
  res.df <- res.df %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    gather(key = "gene",value = "gene.exp",-cell_id,-Time) %>%
    mutate(gene=factor(gene,levels = c(tmp.Lgene,tmp.Rgene)))
  
  
  
  p <- ggplot(data = res.df,aes_string(x="Time",y="gene.exp",colour="gene"))+
    geom_point(alpha=0.1,size=1)+
    geom_smooth(aes_string(group="gene"),se = F)+
    scale_color_manual(values = c("#E41A1C","#377EB8"))+
    ylab("gene.exp")+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 18)+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 18 ) ,
          axis.text = element_text(size = 18),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ),
          axis.text.x = element_text(hjust = 0.8,angle = 30))
  
  return(p)
}



#### Define function
myDavid_GOplot <- function(df=tmp,
                           use.FDR=F,
                           FDR.cut=0.01,
                           p.cut=0.01,
                           show_number=20,
                           title="LR_cor_ge_0.1_Lgene",
                           fill.color="#E64B35FF",
                           plot.ylab="BP",
                           term.pos.adjust=0.8,
                           font.size = 12){
  
  data.plot <- df
  
  data.plot$logp <- -log10(data.plot$PValue)
  data.plot$logFDR <- -log10(data.plot$FDR)
  data.plot$logp.adjust <- -log10(data.plot$Benjamini)
  
  data.plot$Term <- gsub(pattern = ".*~",replacement = "",data.plot$Term)
  
  if(use.FDR){
    data.plot <- data.plot %>%
      filter(FDR < FDR.cut) %>%
      top_n(n = show_number,wt = logFDR)
  }else{
    data.plot <- data.plot %>%
      filter(PValue < p.cut) %>%
      top_n(n = show_number,wt = logp)
  }
  
  data.plot$Term <- factor(data.plot$Term,levels = rev(data.plot$Term))
  data.term <- data.frame(x = 0,
                          y = row_number(data.plot$Term),
                          label = data.plot$Term,
                          stringsAsFactors = F)
  p1 <- ggplot() + 
    geom_bar(data = data.plot, aes_string(x = "Term", 
                                          y = ifelse(use.FDR,"logFDR","logp")),
             stat = "identity",
             fill=fill.color)+
    coord_flip() +
    ggtitle(title)+
    xlab(plot.ylab)+
    ylab(ifelse(use.FDR,"-log10(FDR)","-log10(Pvalue)"))+
    geom_text(data = data.term, 
              aes(x=y,y=x,label = label), 
              color = "black",
              hjust = term.pos.adjust,
              size = 5.8)+
    theme_cowplot(font_size = font.size) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0.5))
  
  return(p1)
}

#### myenrichrGo_plot
myenrichr_GOplot <- function(df=tmp,
                           p.adjust.cut=0.01,
                           show_number=20,
                           title="LR_cor_ge_0.1_Lgene",
                           fill.color="#E64B35FF",
                           plot.ylab="BP",
                           term.pos.adjust=0.8,
                           font.size = 12){
  data.plot <- df
  data.plot$logp.adjust <- -log10(data.plot$p.adjust)
  data.plot$Term <- gsub(pattern = "^GO_",replacement = "",data.plot$ID)
  data.plot <- data.plot %>%
    dplyr::filter(p.adjust < p.adjust.cut) %>%
    slice_max(order_by = logp.adjust,n=show_number,with_ties=F)
    
    
  data.plot$Term <- factor(data.plot$Term,levels = rev(data.plot$Term))
  data.term <- data.frame(x = 0,
                          y = row_number(data.plot$Term),
                          label = data.plot$Term,
                          stringsAsFactors = F)
  p1 <- ggplot() + 
    geom_bar(data = data.plot, aes_string(x = "Term", 
                                          y = "logp.adjust"),
             stat = "identity",
             fill=fill.color)+
    coord_flip() +
    ggtitle(title)+
    xlab(plot.ylab)+
    ylab("-log10(p.adjust)")+
    geom_text(data = data.term, 
              aes(x=y,y=x,label = label), 
              color = "black",
              hjust = term.pos.adjust,
              size = 5.8)+
    theme_cowplot(font_size = font.size) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0.5))
  
  return(p1)
}





myCreateSeurat <- function (data=data, project = "SeuratProject", assay = "RNA", min.cells = 0, 
                            min.features = 0, names.field = 1, names.delim = "_", meta.data = NULL) {
  if (!is.null(x = meta.data)) {
    if (is.null(x = rownames(x = meta.data))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(x = rownames(x = meta.data), y = colnames(x = data)))) {
      warning("Some cells in meta.data not present in provided data matrix.")
      meta.data <- meta.data[intersect(x = rownames(x = meta.data), 
                                       y = colnames(x = data)), ]
    }
    if (is.data.frame(x = meta.data)) {
      new.meta.data <- data.frame(row.names = colnames(x = data))
      for (ii in 1:ncol(x = meta.data)) {
        new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <- meta.data[, 
                                                                                         ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  }
  assay.data <- CreateAssayObject(data=as.matrix(data))
  Key(object = assay.data) <- paste0(tolower(x = assay), "_")
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay
  init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay]]))
  idents <- factor(x = unlist(x = lapply(X = colnames(x = assay.data), 
                                         FUN = Seurat:::ExtractField, field = names.field, delim = names.delim)))
  if (any(is.na(x = idents))) {
    warning("Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name")
  }
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == 
      length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = assay.data))
  }
  names(x = idents) <- colnames(x = assay.data)
  object <- new(Class = "Seurat", assays = assay.list, meta.data = init.meta.data, 
                active.assay = assay, active.ident = idents, project.name = project, 
                version = packageVersion(pkg = "Seurat"))
  object[["orig.ident"]] <- idents
  n.calc <- Seurat:::CalcN(object = assay.data)
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, 
                               sep = "_")
    object[[names(x = n.calc)]] <- n.calc
  }
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

myRemoveNA <- function(df){
  df[is.na(df)] <- 0
  return(df)
}


myNormalize <- function(data.exp){
  data.exp <- log2(data.exp+1)
  return(data.exp)
}

myMatrixFilter <- function(data.exp,filter=0.95){
  data.exp <- apply(data.exp,2,FUN = function(x) ifelse(x >= quantile(x,filter),quantile(x,filter),x))
  return(data.exp)
}

myBoxPlot <- function(data.plot,x,y,xlab,ylab,...){
  p <- ggplot(data.plot,aes_string(x,y,...))+
    geom_boxplot()+
    theme_cowplot(font_size = 15)+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.position = "right",
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ),
          axis.text.x = element_text(hjust = 0.8,
                                     angle = 30))
  return(p)
}

### file generate

myFileName <- function(prefix,suffix){
  res <- paste0(prefix,"_",format(Sys.Date(),"%Y%m%d"),suffix)
  return(res)
}

### Capitlized

myCapitalized <- function(tmp.string){
  tmp.res <- Hmisc::capitalize(tolower(tmp.string))
}



### default to save to jpg
myggsave <- function(p,prefix,suffix,width=8,height=8,...){
  ggsave(p,filename = myFileName(prefix = prefix,suffix = suffix),width = width,height = height,...)
}


myrunDAVID <- function (gene=eg$ENTREZID, 
                        idType = "ENTREZ_GENE_ID", 
                        annotation = "GOTERM_BP_DIRECT",
                        david.user="longtengwang@pku.edu.cn") {
  
  david.pkg <- "RDAVIDWebService"
  pkgs <- installed.packages()[, 1]
  if (!david.pkg %in% pkgs) {
    stop("You should have RDAVIDWebService package installed before using enrichDAVID...")
  }
  require(david.pkg, character.only = TRUE)
  
  david <- DAVIDWebService$new(email = david.user, 
                               url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  idType <- match.arg(idType, getIdTypes(david))
  david.res <- addList(david, gene, idType = idType, listName = "clusterProfiler", 
                       listType = "Gene")
  if (david.res$inDavid == 0) {
    stop("All id can not be mapped. Please check 'idType' parameter...")
  }
  
  setAnnotationCategories(david, annotation)
  x <- getFunctionalAnnotationChart(david)
  if (length(x@.Data) == 0) {
    warning("No significant enrichment found...")
    return(NULL)
  }
  term.df <- x@.Data
  names(term.df) <- x@names
  term.df <- as.data.frame(term.df)
  return(term.df)
}

myPerformDavidGO <- function(gene.lits,output.table,output.fig,plot.title,annotation){
  eg = bitr(gene.lits, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  term.df = myrunDAVID(gene = eg$ENTREZID, 
                       idType="ENTREZ_GENE_ID",
                       annotation=annotation,
                       david.user = "longtengwang@pku.edu.cn")
  term.df$Genes <- as.character(term.df$Genes)
  
  test <- lapply(term.df$Genes,function(x){
    eg <- bitr(unlist(strsplit(x, split = ", ")),
               fromType = "ENTREZID",
               toType = "SYMBOL",
               OrgDb = org.Mm.eg.db)
    res <- paste(eg$SYMBOL,collapse = ", ")
    return(res)
  })
  term.df$Genes <- unlist(test)
  write.table(term.df,file = output.table,sep = "\t",row.names = F,quote = F)
  
  p <- myDavid_GOplot(df = term.df,
                      use.FDR = F,
                      p.cut = 0.01,
                      plot.ylab = annotation,
                      title = plot.title,
                      show_number = 20)
  ggsave(p,filename = output.fig,width = 8,height = 8,dpi = 300)
}


myMatToZeroOne <- function(mat,cutoff){
  mat <- ifelse(mat > cutoff,1,0)
  return(mat)
}

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

myMatCutoff <- function(df,cutoff){
  res <- df[apply(df,1,function(x) all(x>cutoff)),]
  return(res)
}

median_center <- function(x){
  res <- x - apply(x, 2, median)
  return(res)
}


myPerformDavidGO <- function(gene.lits,output.table,output.fig,plot.title,annotation){
  eg = bitr(gene.lits, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  term.df = myrunDAVID(gene = eg$ENTREZID, 
                       idType="ENTREZ_GENE_ID",
                       annotation=annotation,
                       david.user = "longtengwang@pku.edu.cn")
  term.df$Genes <- as.character(term.df$Genes)
  
  test <- lapply(term.df$Genes,function(x){
    eg <- bitr(unlist(strsplit(x, split = ", ")),
               fromType = "ENTREZID",
               toType = "SYMBOL",
               OrgDb = org.Mm.eg.db)
    res <- paste(eg$SYMBOL,collapse = ", ")
    return(res)
  })
  term.df$Genes <- unlist(test)
  write.table(term.df,file = output.table,sep = "\t",row.names = F,quote = F)
  
  p <- myDavid_GOplot(df = term.df,
                      use.FDR = F,
                      p.cut = 0.01,
                      plot.ylab = annotation,
                      title = plot.title,
                      show_number = 20)
  ggsave(p,filename = output.fig,width = 8,height = 8,dpi = 300)
}

myBoxplot_advanced <- function(data = df,
                               x = "dose",
                               y = "len",
                               fill = "supp",
                               color = "supp",
                               mycolor = mypalette,
                               boxwidth = 0.6,
                               errorbarwidth = 0.25,
                               dodgingwidth = 0.8,
                               ylimts = c(0,40),
                               ybreaks = seq(0,40,10),
                               fontsize = 18,
                               linesize = 1,
                               medianlinesize = 1,
                               axislinesize = 1,
                               outlier.shape = 19,
                               outlier.alpha = 0.1){
  bxp <- ggplot(data = data,aes_string(x = x,
                                       y = y,
                                       fill = fill,
                                       color= color))+
    geom_boxplot(linetype = "dashed",
                 position = position_dodge(dodgingwidth),
                 width = boxwidth, 
                 size = linesize,
                 outlier.shape = outlier.shape,
                 outlier.alpha = outlier.alpha) +
    stat_boxplot(geom = "errorbar", 
                 aes(ymin = ..ymax..),
                 width = errorbarwidth , 
                 size = linesize,
                 position = position_dodge(dodgingwidth)) +
    stat_boxplot(geom = "errorbar", 
                 aes(ymax = ..ymin..),
                 width = errorbarwidth ,
                 position = position_dodge(dodgingwidth)) +
    stat_boxplot(aes(ymin = ..lower.., 
                     ymax = ..upper..),
                 outlier.alpha = outlier.alpha,
                 outlier.shape = outlier.shape,
                 width = boxwidth,
                 size = linesize,
                 position = position_dodge(dodgingwidth)) +
    ### The bang-bang operator !! forces a single object. One common case for  !! is to substitute an environment-variable (created with <-) with a data-variable (inside a data frame).
    stat_boxplot(data = data,geom = "errorbar", 
                 aes(x = !!sym(x),
                     y = !!sym(y),
                     ymin = ..middle..,
                     ymax = ..middle..),
                 size = medianlinesize,
                 col="black",
                 width = boxwidth ,
                 position = position_dodge(dodgingwidth))
  if(is.null(ylimts)){
    bxp <- bxp + scale_fill_manual(values = mycolor)+
      scale_color_manual(values =mycolor)+
      theme_cowplot(font_size = fontsize)+
      theme(legend.position = "right",
            axis.line.x = element_line( size = axislinesize ),
            axis.line.y = element_line( size = axislinesize ),
            axis.text.x = element_text(hjust = 0.8,
                                       angle = 30))
  }else{
    bxp <- bxp + scale_y_continuous(limits = ylimts,breaks = ybreaks,expand = expansion(mult = c(0,0.1)))+
      scale_fill_manual(values = mycolor)+
      scale_color_manual(values = mycolor)+
      theme_cowplot(font_size = fontsize)+
      theme(legend.position = "right",
            axis.line.x = element_line( size = axislinesize ),
            axis.line.y = element_line( size = axislinesize ),
            axis.text.x = element_text(hjust = 0.8,
                                       angle = 30))
  }
  return(bxp)
}

#### about color bar

# code from
# https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
# edited
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab=NA, yaxt='n', ylab=NA, main=title)
  #axis(side = 4, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# code from
# https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r 
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#### myshowcolorbar

my_showcolorbar <- function(tmp.colorbar,tmp.name,n=100){
  image(
    1:n, 1, as.matrix(1:n),
    col = tmp.colorbar(n),
    xlab = tmp.name, 
    ylab = "", 
    xaxt = "n", 
    yaxt = "n", 
    bty = "n"
  )
}

#####started to edit this Utils function
#####202108
#####20210927
#####20211023

#####---------1. Operation function--------

myGetGRangesFromsTxDb <- function(
  txdb = mm9KG_txdb,
  standard.chromosomes = TRUE,
  verbose = TRUE
) {
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n",
         "https://www.bioconductor.org/packages/biovizBase/")
  }
  # convert seqinfo to granges
  whole.genome <-  as(object = seqinfo(x = txdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
  }
  # extract genes from each chromosome
  if (verbose) {
    ####crunch featching GRanges from various data source
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      cat(seqlevels(whole.genome)[x],sep = "\n")
      biovizBase::crunch(
        obj = txdb,
        which = whole.genome[x],
        columns = c("tx_id", "tx_name", "gene_id"))
    })
  } else {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      suppressMessages(expr = biovizBase::crunch(
        obj = txdb,
        which = whole.genome[x],
        columns = c("tx_id", "tx_name","gene_id")))
    })
  }
  # combine
  tx <- do.call(what = c, args = tx)
  return(tx)
}


p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}

my_enrichment.pvalue <- function(M,N,k,n,lower.tail=F){
  p <- phyper(k-1,M,N-M,n,lower.tail = lower.tail)
  return(p)
}
myTrunclog <- function(x){
  res <- ifelse(x>1,log10(x),0)
}

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


#####---------2. IO--------

myBegraphToBigwig <- function(tmp.input,tmp.output,region=NULL,tmp.seqlength=NULL,tmp.levels=NULL){
  cat("load bedgraph",sep = "\n")
  if(!is.null(region)){
    region_data <- rtracklayer::import.bedGraph(
      con = tmp.input,
      which = as(region,"GRanges")
    )
  }else{
    region_data <- rtracklayer::import.bedGraph(
      con = tmp.input,
    )
  }
  #rtracklayer::close(con = tmp.input)
  seqlevels(region_data) <- tmp.levels
  seqlengths(region_data) <- tmp.seqlength
  cat("export bedgraph",sep = "\n")
  rtracklayer::export.bw(object = region_data,con = tmp.output)
  #rtracklayer::close(con = tmp.output)
}

#####---------3. visualization--------
mycolor.bar <- function(my.colors,min,max=-min,vertical=T,showticks=T,tmp.title=NULL){
  
  if(vertical){
    z=matrix(1:100,nrow=1)
    x=1
    y=seq(min,max,len=100) # supposing 3 and 2345 are the range of your data
    image(x,y,z,col=my.colors,axes=FALSE,xlab="",ylab="")
    axis(4,tick = showticks,lwd.ticks = 1)
    title(main = tmp.title)
  }
  else{
    x = seq(min,max,len=100)
    y = 1
    z=matrix(1:100,ncol=1)
    image(x,y,z,col=my.colors,axes=FALSE,xlab="",ylab="")
    axis(1,tick = showticks,lwd.ticks = 1)
    title(main = tmp.title)
  }
}


myBigwigTrack <- function(
  region = tmp.region,
  bigwig = K4me3,
  smooth = 200,
  lognorm = F,
  type = "coverage",
  y_label = "Score",
  fontsize=18,
  track.color="blue",
  tmp.ylimits=c(0,16),
  max.downsample = 3000,
  downsample.rate = 0.1,
  tmp.seed=42
) {
  possible_types <- c("line", "heatmap", "coverage","bar")
  if (!(type %in% possible_types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible_types, collapse = ", ")
    )
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- as(region,"GRanges")
  }
  ### to smooth track data
  region_data <- rtracklayer::import(
    con = bigwig,
    which = region,
    as = "NumericList"
  )[[1]]
  if (!is.null(x = smooth)) {
    region_data <- RcppRoll::roll_mean(x = region_data, n = smooth, fill = 0L)
  }
  region_data <- data.frame(
    position = start(x = region):end(x = region),
    score = region_data,
    stringsAsFactors = FALSE
  )
  window.size = width(x = region)
  sampling <- min(max.downsample, window.size * downsample.rate)
  set.seed(tmp.seed)
  coverages <- slice_sample(.data = region_data, n = sampling)
  if(lognorm){
    coverages$score <- myTrunclog(coverages$score)
  }
  p <- ggplot(
    data = coverages,
    mapping = aes_string(x = "position", y = "score")
  )
  if (type == "line") {
    p <- p + geom_line()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    region_data$bin <- floor(x = region_data$position / smooth)
    region_data <- group_by(region_data, bin)
    region_data <- mutate(region_data, score = mean(x = score))
    region_data <- ungroup(region_data)
    region_data <- unique(x = region_data[, c("bin", "score")])
    p <- ggplot(
      data = region_data,
      mapping = aes_string(x = "bin", y = 1, fill = "score")
    ) + geom_tile() + scale_fill_viridis_c()
  } else if (type == "coverage") {
    p <- p + 
      geom_area(fill=track.color,color=track.color)
  } else if (type == "bar") {
    p <- p + 
      geom_bar(stat="identity",fill=track.color,color=track.color)
  } 
  chromosome <- as.character(x = seqnames(x = region))
  p <- p  +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label) +
    scale_x_continuous(expand = c(0,0),limits = c(start(region),end(region)))+
    scale_y_continuous(expand = c(0,0),limits = tmp.ylimits,
                       breaks = c(tmp.ylimits[2]),
                       labels = paste0("[",0,"-",tmp.ylimits[2],"]"))+
    theme_cowplot(font_size = fontsize)+
    theme(axis.line = element_line(size = 0.5),
          axis.ticks = element_line(size = 0.5),
          axis.title.y = element_text(angle = 0),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
  return(p)
}


split_body <- function(df, width = 1000) {
  wd <- df$end - df$start
  nbreak <- wd / width
  if (nbreak > 1) {
    steps <- 0:(nbreak)
    starts <- (width * steps) + df$start
    starts[starts > df$end] <- NULL
  } else {
    starts <- df$end
  }
  breaks <- data.frame(
    seqnames = df$seqnames[[1]],
    start = starts,
    end = starts + 1,
    strand = df$strand[[1]],
    gene_name = df$gene_name[[1]],
    type = "arrow"
  )
  return(breaks)
}

record_overlapping <- function(annotation = gene_bodies, min.gapwidth = 1000) {
  # convert back to granges
  annotation.stash <- annotation
  annotation$strand <- "*"
  gr <- makeGRangesFromDataFrame(
    df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
  )
  # work out which ranges overlap
  collapsed <- GenomicRanges::reduce(
    x = gr, with.revmap = TRUE, min.gapwidth = min.gapwidth
  )$revmap
  idx <- seq_along(gr)
  for (i in seq_along(collapsed)) {
    mrg <- collapsed[[i]]
    for (j in seq_along(mrg)) {
      idx[[mrg[[j]]]] <- j
    }
  }
  names(x = idx) <- gr$gene_name
  return(idx)
}

reformat_annotations <- function(
  annotation = annotation.subset,
  start.pos = start.pos,
  end.pos = end.pos,
  arrow_sbreaks=2000
) {
  annotation <- annotation[annotation$type == "exon"]
  exons <- as.data.frame(x = annotation)
  annotation <- split(
    x = annotation,
    f = as.character(annotation$gene_id)
  )
  annotation <- lapply(X = annotation, FUN = as.data.frame)
  
  # add gene total start / end
  gene_bodies <- lapply(seq_along(annotation), function(i){
    df <- data.frame(
      seqnames = annotation[[i]]$seqnames[[1]],
      start = min(annotation[[i]]$start),
      end = max(annotation[[i]]$end),
      strand = annotation[[i]]$strand[[1]],
      gene_name = annotation[[i]]$gene_id[[1]],
      type = "body"
    )
    # trim any that extend beyond region
    df$start <- ifelse(
      test = df$start < start.pos,
      yes = start.pos,
      no = df$start
    )
    df$end <- ifelse(
      test = df$end > end.pos,
      yes = end.pos,
      no = df$end
    )
    breaks <- split_body(df = df,width = arrow_sbreaks)
    df <- rbind(df, breaks)
    return(df)
  })
  
  gene_bodies <- do.call(what = rbind, args = gene_bodies)
  
  # record if genes overlap
  overlap_idx <- record_overlapping(annotation = gene_bodies, min.gapwidth = 1000)
  gene_bodies$dodge <- overlap_idx[as.character(gene_bodies$gene_name)]
  exons$dodge <- overlap_idx[as.character(exons$gene_id)]
  
  label_df <- gene_bodies[gene_bodies$type == "body", ]
  label_df$width <- label_df$end - label_df$start
  label_df$position <- label_df$start + (label_df$width / 2)
  
  onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
  onminus <- gene_bodies[gene_bodies$strand == "-", ]
  
  return(
    list(
      "labels" = label_df,
      "exons" = exons,
      "plus" = onplus,
      "minus" = onminus
    )
  )
}



myGenePlot <- function(annotation=tx, 
                       region="chr7:152042291-152056148",showlabel=T,
                       arrow_sbreaks=2000,font_size=18,label_size=8) {
  
  if (is.null(x = annotation)) {
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- as(region,"GRanges")
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  genes.keep <- unique(x = annotation.subset$gene_id)
  annotation.subset <- annotation[
    fastmatch::fmatch(x = annotation$gene_id, table = genes.keep, nomatch = 0L) > 0L
  ]
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      arrow_sbreaks=arrow_sbreaks
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 5
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.05, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.05, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + (n_stack * 0.2)
    if(showlabel){
      p <- p + geom_text(
        data = annotation_df_list$labels,
        mapping = aes_string(x = "position", y = "dodge", label = "gene_name"),
        size = label_size,fontface="bold"
      )
    }
    y_limit <- c(0.9, n_stack + (n_stack * 0.5))
  }
  p <- p +
    theme_cowplot(font_size = font_size) +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    # scale_x_continuous(expand = c(0,0),
    #                    limits = c(start.pos,end.pos),
    #                    breaks = seq(start.pos,end.pos,length.out=5))+
    scale_x_continuous(expand = c(0,0),
                       limits = c(start.pos,end.pos))+
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}


get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

#####---------3.1 defien colors --------
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




