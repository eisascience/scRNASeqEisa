# Customized functions from or for Seurat analysis...


AddModuleScoreSindv <- function(Ser_object=NULL, geneSet=NULL){
  # compute individual scores of genes in a gene list 
  #      as oppose to Seurat default ... it too does it individually, but controls the naming
  
  # Also we scale the results between 0 and 1 for probability equivalence 
  #     the range01() is in the Useful.R function set
  
  
  for(geneN in geneSet){
    print(geneN)
    #geneN = geneSet[1]
    Ser_object <- AddModuleScore(Ser_object, 
                                 genes.list = geneN, 
                                 genes.pool = rownames(Ser_object@data), 
                                 n.bin = 25,
                                 seed.use = 1, 
                                 ctrl.size = 100, 
                                 use.k = FALSE, 
                                 enrich.name = geneN,
                                 random.seed = 1)
    
    Ser_object@meta.data[,paste(geneN, "1", sep="")] <- range01(Ser_object@meta.data[,paste(geneN, "1", sep="")])
    
  }
  #colnames(Ser_object@meta.data[,paste(geneSet, "1", sep="")]) <- geneSet
  return(Ser_object)
} 


AddModuleScoreV2 <- 
  function (object, genes.list = NULL, genes.pool = NULL, n.bin = 25, 
            seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name = "Cluster", 
            random.seed = 1, returnSerObj=T) 
  {
    # object = CleaningLS$SeurComboObj
    # random.seed = 1
    # n.bin = 25
    # genes.list = genes.list <- list(S.Score = s.genes, G2M.Score = g2m.genes)
    # enrich.name <- "Cell Cycle"
    # ctrl.size = min(vapply(X = genes.list, 
    #            FUN = length, FUN.VALUE = numeric(1)))
    # genes.pool = NULL
    
    
    set.seed(seed = random.seed)
    genes.old <- genes.list
    # use.k = FALSE
    
    # if (use.k) {
    #   genes.list <- list()
    #   for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
    #     genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
    #                                          i))
    #   }
    #   cluster.length <- length(x = genes.list)
    # }
    # else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(X = genes.list, FUN = function(x) {
      return(intersect(x = x, y = rownames(x = object@data)))
    })
    
    cluster.length <- length(x = genes.list)
    #}
    # if (!all(LengthCheck(values = genes.list))) {
    #   warning(paste("Could not find enough genes in the object from the following gene lists:", 
    #                 paste(names(x = which(x = !LengthCheck(values = genes.list)))), 
    #                 "Attempting to match case..."))
    #   genes.list <- lapply(X = genes.old, FUN = CaseMatch, 
    #                        match = rownames(x = object@data))
    # }
    # if (!all(LengthCheck(values = genes.list))) {
    #   stop(paste("The following gene lists do not have enough genes present in the object:", 
    #              paste(names(x = which(x = !LengthCheck(values = genes.list)))), 
    #              "exiting..."))
    # }
    
    if (is.null(x = genes.pool)) {
      genes.pool = rownames(x = object@data)
    }
    
    data.avg <- Matrix::rowMeans(x = object@data[genes.pool, ])
    data.avg <- data.avg[order(data.avg)]
    data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
    names(x = data.cut) <- names(x = data.avg)
    ctrl.use <- vector(mode = "list", length = cluster.length)
    for (i in 1:cluster.length) {
      genes.use <- genes.list[[i]]
      for (j in 1:length(x = genes.use)) {
        ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                                data.cut[genes.use[j]])], size = ctrl.size, replace = FALSE)))
      }
    }
    ctrl.use <- lapply(X = ctrl.use, FUN = unique)
    ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                          ncol = ncol(x = object@data))
    for (i in 1:length(ctrl.use)) {
      genes.use <- ctrl.use[[i]]
      ctrl.scores[i, ] <- Matrix::colMeans(x = object@data[genes.use, 
                                                           ])
    }
    genes.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                           ncol = ncol(x = object@data))
    for (i in 1:cluster.length) {
      genes.use <- genes.list[[i]]
      data.use <- object@data[genes.use, , drop = FALSE]
      genes.scores[i, ] <- Matrix::colMeans(x = data.use)
    }
    genes.scores.use <- genes.scores - ctrl.scores
    rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
    genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
    rownames(x = genes.scores.use) <- colnames(x = object@data)
    
    if(returnSerObj){
      object <- AddMetaData(object = object, metadata = genes.scores.use, 
                            col.name = colnames(x = genes.scores.use))
      gc(verbose = FALSE)
      return(object)
    } else {
      return(genes.scores.use)
      
    }
    
  }



ggUMAP <- function(object, colFac = NULL, col_vector=NULL, ptSize=0.15, ptAlpha=0.5, add2title="", 
                   legendTitle="", cells.use = NULL){
  #updated Feb/19/19 - EM :: use.cell is needed to subset cells
  
  
  datat.temp <- as.data.frame(object@dr$umap@cell.embeddings)
  if(!is.null(cells.use)){
    datat.temp <- datat.temp[cells.use, ]
  }
  
  if(!is.null(colFac)){
    
    if(!is.factor(colFac)) colFac <- factor(colFac)
    
    if(is.null(col_vector)) col_vector = gg_color_hue(length(levels(colFac)))
    
    
    myColors <- col_vector[1:length(levels(colFac))]
    names(myColors) <- levels(colFac)
    
    colScale <- scale_colour_manual(name = ifelse(legendTitle=="", "grp", legendTitle), values = myColors)
    
    
    datat.temp$colFac <- colFac
    
    ggplot(datat.temp, aes(UMAP1, UMAP2 , col=colFac)) + 
      geom_point(size=ptSize, alpha = ptAlpha) +
      theme(legend.position = "bottom") +
      ggtitle(paste("Expression of factor on UMAP", add2title, sep="")) + 
      theme_bw() + colScale + guides(colour = guide_legend(override.aes = list(size=8, alpha=1)))
    
  } else {
    
    ggplot(datat.temp, aes(UMAP1, UMAP2)) + 
      geom_point(size=ptSize, alpha = ptAlpha) +
      theme(legend.position = "bottom") +
      ggtitle(paste("UMAP 2D Space", add2title, sep="")) + 
      theme_bw() 
    
  }
  
}




plotGeneExprUMAP <- function(datat, DGEmat, geneName, logExprBL = T, ThresholdExpr = F, ExtraTitle=""){
  #updated Feb/21/19 - EM :: ExtraTitle and theme adjustment
  
  library(viridis)
  datat.temp <- datat
  datat.temp$GeneExpr <- DGEmat[geneName,rownames(datat.temp)]
  
  if(ThresholdExpr!=F){
    datat.temp$GeneExpr <- ifelse(asinh(datat.temp$GeneExpr)>=ThresholdExpr, 1, 0)
    
  }
  if(logExprBL) ggp1 <- ggplot(datat.temp, aes(UMAP1, UMAP2 , col=log(GeneExpr+1)))
  if(!logExprBL) ggp1 <- ggplot(datat.temp, aes(UMAP1, UMAP2 , col=asinh(GeneExpr)))
  ggp1 + geom_point(size=0.15, alpha = .5) +
    scale_color_viridis(direction=1) +
    theme(legend.position = "bottom") +
    ggtitle(paste("Expression of gene ", geneName,  " on UMAP\n",ExtraTitle, sep="")) + theme_bw()
}

plotMultiGeneExprUMAP <- function(datat, DGEmat, geneName, logExprBL = T, ThresholdExpr = F, ExtraTitle=""){
  #updated Feb/21/19 - EM :: ExtraTitle and theme adjustment
  
  library(viridis)
  datat.temp <- datat
  if(length(geneName)==1)   datat.temp$GeneExpr <- DGEmat[geneName,rownames(datat.temp)]
  if(length(geneName)>1)    datat.temp$GeneExpr <- Matrix::colSums(DGEmat[geneName,rownames(datat.temp)])
  
  if(ThresholdExpr!=F){
    datat.temp$GeneExpr <- ifelse(asinh(datat.temp$GeneExpr)>=ThresholdExpr, 1, 0)
    
  }
  if(logExprBL) ggp1 <- ggplot(datat.temp, aes(UMAP1, UMAP2 , col=log(GeneExpr+1)))
  if(!logExprBL) ggp1 <- ggplot(datat.temp, aes(UMAP1, UMAP2 , col=asinh(GeneExpr)))
  ggp1 + geom_point(size=0.15, alpha = .5) +
    scale_color_viridis(direction=1) +
    theme(legend.position = "bottom") +
    ggtitle(paste("Expression of sum of ", length(geneName),  "genes on UMAP\n",ExtraTitle, sep="")) + theme_bw() 
}


FindMarkers_CellType <- function(Ser_object, celltype, celltype2 = NULL, assay.type = "RNA", TriTestMode=T){
  #Updated Feb/25/2019 celltype2
  
  cells.1 <- WhichCells(object = Ser_object, ident = celltype)
  if(is.null(celltype2)){
    cells.2 <- WhichCells(object = Ser_object, ident.remove = celltype)
    cells.2 <- setdiff(x = cells.2, y = cells.1)
  } else {
    cells.2 <- WhichCells(object = Ser_object, ident = celltype2)
  }
  
  
  data.use <- GetAssayData(object = Ser_object, 
                           assay.type = assay.type, 
                           slot = "data")
  genes.use <- rownames(x = data.use)
  
  
  data.temp1 <- round(x = apply(X = data.use[genes.use, cells.1, 
                                             drop = F], MARGIN = 1, FUN = function(x) {
                                               return(sum(x > 0)/length(x = x))
                                             }), digits = 3)
  
  data.temp2 <- round(x = apply(X = data.use[genes.use, cells.2, 
                                             drop = F], MARGIN = 1, FUN = function(x) {
                                               return(sum(x > 0)/length(x = x))
                                             }), digits = 3)
  
  
  data.alpha <- cbind(data.temp1, data.temp2)
  colnames(x = data.alpha) <- c("pct.1", "pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  genes.use <- names(x = which(x = alpha.min > 0.1))
  
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, 
                                  FUN = min)
  genes.use <- names(x = which(x = alpha.min > .1 & alpha.diff > 
                                 -Inf))
  
  
  data.1 <- apply(X = data.use[genes.use, cells.1, drop = F], 
                  MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 
                                                      1))
  data.2 <- apply(X = data.use[genes.use, cells.2, drop = F], 
                  MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 
                                                      1))
  total.diff <- (data.1 - data.2)
  
  genes.diff <- names(x = which(x = abs(x = total.diff) > 
                                  0.25))
  
  genes.use <- intersect(x = genes.use, y = genes.diff)
  
  to.return <- WilcoxDETest(object = Ser_object, 
                            cells.1 = cells.1, 
                            cells.2 = cells.2, 
                            genes.use = genes.use, 
                            print.bar = T, 
                            assay.type=assay.type)
  colnames(to.return) <- c("p_val_Wilcox")
  
  if(TriTestMode){
    to.return_bimod <- DiffExpTest(object = Ser_object, 
                                   assay.type = assay.type, 
                                   cells.1 = cells.1, 
                                   cells.2 = cells.2, 
                                   genes.use = genes.use, 
                                   print.bar = T)
    colnames(to.return_bimod) <- c("p_val_bimod")
    
    to.return_roc <- MarkerTest(object = Ser_object, 
                                assay.type = assay.type, 
                                cells.1 = cells.1, 
                                cells.2 = cells.2, 
                                genes.use = genes.use, 
                                print.bar = T)
    colnames(to.return_roc) <- c("p_val_roc")
  }
  
  
  
  to.return[, "avg_logFC"] <- total.diff[rownames(x = to.return)]
  
  if(TriTestMode) to.return <- cbind(to.return, to.return_bimod, to.return_roc, data.alpha[rownames(x = to.return), , drop = FALSE])
  if(!TriTestMode) to.return <- cbind(to.return, data.alpha[rownames(x = to.return), , drop = FALSE])
  
  to.return$p_val_Wilcox_adj_bonf = p.adjust(p = to.return$p_val_Wilcox, method = "bonferroni", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                                          slot = "data")))
  
  to.return$p_val_Wilcox_adj_fdr = p.adjust(p = to.return$p_val_Wilcox, method = "fdr", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                                  slot = "data")))
  
  to.return$p_val_Wilcox_adj_BH = p.adjust(p = to.return$p_val_Wilcox, method = "BH", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                                slot = "data")))
  
  if(TriTestMode) {
    to.return$p_val_bimod_adj_bonf = p.adjust(p = to.return$p_val_bimod, method = "bonferroni", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                                          slot = "data")))
    
    to.return$p_val_bimod_adj_fdr = p.adjust(p = to.return$p_val_bimod, method = "fdr", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                                  slot = "data")))
    
    to.return$p_val_bimod_adj_BH = p.adjust(p = to.return$p_val_bimod, method = "BH", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                                slot = "data")))
    
    
    to.return$p_val_roc_adj_bonf = p.adjust(p = to.return$p_val_roc, method = "bonferroni", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                                      slot = "data")))
    
    to.return$p_val_roc_adj_fdr = p.adjust(p = to.return$p_val_roc, method = "fdr", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                              slot = "data")))
    
    to.return$p_val_roc_adj_BH = p.adjust(p = to.return$p_val_roc, method = "BH", n = nrow(x = GetAssayData(object = Ser_object, assay.type = assay.type, 
                                                                                                            slot = "data")))
    
  }
  
  to.return <- as.data.frame(to.return[order(to.return$p_val_Wilcox, -to.return$avg_logFC),])
  
  return(to.return)
}





