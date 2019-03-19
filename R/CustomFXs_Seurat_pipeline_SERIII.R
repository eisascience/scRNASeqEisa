hasStepRun <- function(seuratObj, name) {
  return(!is.null(seuratObj@misc[[paste0(name, 'Run')]]))
}

markStepRun <- function(seuratObj, name, saveFile = NULL) {
  seuratObj@misc[paste0(name, 'Run')] <- T
  if (!is.null(saveFile)){
    saveRDS(seuratObj, file = saveFile)
  }

  return(seuratObj)
}


createSeuratObj_SERIII <- function(seuratData = NA, project = NA, minFeatures = 0, minCells=0){
  #BBimber
  seuratObj <- CreateSeuratObject(counts = seuratData,
                                  min.cells = minCells,
                                  min.features = minFeatures,
                                  project = project)

  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObj), value = TRUE)
  p.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
  seuratObj[['p.mito']] <- p.mito


  return(seuratObj)
}


MakeSerObjs_10XFolders_SERIII <- function(counts.path = NULL,
                                   min.cells = 0,
                                   min.genes = 0,
                                   ProjName="10X",
                                   save.path = NULL,
                                   returnList=F, path.exclude="raw", string.exclude=NULL){

  #this function searches for 10X folder.
  #this is an update from the MakeSerObjs_10XFolders() which was done in Seurat II
  #unfortuantely, as much as the developers know what a pain it is to change syntax
  # when going from one version to anther, they STILL chose to do this :(


  require(Seurat)
  if(returnList) TempLS <- list()

  if(!is.null(counts.path)){
    if(is.null(save.path)) save.path <- counts.path

    exp.dirs <- list.files(counts.path, recursive = T, full.names = T, pattern = ".mtx")
    exp.dirs <- exp.dirs[!grepl(".gz", exp.dirs)]
    exp.dirs <- gsub("/matrix.mtx", "", exp.dirs)
    if(!path.exclude=="") exp.dirs <- exp.dirs[!grepl(path.exclude, exp.dirs)]

    FileNames2Save <- exp.dirs
    if(is.null(string.exclude)) string.exclude <- paste(counts.path, "/", sep="")

    if(!is.null(string.exclude)){
      string.exclude <- c(string.exclude, paste(counts.path, "/", sep=""))
      for(iN in 1:length(string.exclude)){
        FileNames2Save <- gsub(string.exclude[iN], "", FileNames2Save)
      }
    }

    print("Found files... here some..")
    head(exp.dirs)

    if(returnList) TempLS$exp.dirs <- exp.dirs

    for(xN in 1:length(exp.dirs)){
      # xN = 1
      # use the list.files below to make sure the expected files are there....
      #list.files(exp.dirs[xN], full.names = T, recursive = T)

      if(returnList) TempLS$SeuratObjs <- list()

      print(exp.dirs[xN])
      print(paste(save.path,"/",FileNames2Save[xN], "_SeuratObj.rds", sep=""))

      if(!file.exists(paste(save.path,"/",FileNames2Save[xN], "_SeuratObj.rds", sep=""))){

        print("Reading in 10X folder...")
        Seurat10X  <- Read10X(data.dir = exp.dirs[xN])

        print("Converting to Seurat Obj....")

        SeuratObjs <- createSeuratObj_SERIII(seuratData = Seurat10X,
                                             minCells = min.cells,   #genes expressed in >= 5 cells
                                             minFeatures = min.genes, #Keep all cells with at least 200 detected genes
                                             project = paste(ProjName, FileNames2Save[xN], sep="_"))



        if(returnList) TempLS$SeuratObjs[[basename(exp.dirs[xN])]] <- SeuratObjs

        print("saving...")
        saveRDS(SeuratObjs,
                paste(save.path,"/",FileNames2Save[xN], "_SeuratObj.rds", sep=""))

      } else {
        print("already Exists...")
      }



    }; remove(xN)

    if(returnList) return(TempLS)
  }

}



printQcPlots1_SERIII <- function(seuratObj) {
  print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))

  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))

  #10x-like plot
  nUMI <- Matrix::colSums(GetAssayData(object = seuratObj, slot = "counts"))
  nUMI <- sort(nUMI)

  countAbove <-unlist(lapply(nUMI, function(x){
    sum(nUMI >= x)
  }))

  plot(log(countAbove), log(nUMI), pch=20, ylab = "UMI/Cell", xlab = "# Cells")
}


processSeurat1_SERIII <- function(seuratObj, mean.cutoff = c(0.0125, 3),
                                  dispersion.cutoff = c(0.5, Inf)){
  #seuratObj = SeurObj_RDS

  if (!hasStepRun(seuratObj, 'NormalizeData')) {
    seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize")
    seuratObj <- markStepRun(seuratObj, 'NormalizeData', saveFile)
  }


  if (!hasStepRun(seuratObj, 'FindVariableFeatures')) {
    seuratObj <- FindVariableFeatures(object = seuratObj,
                                      selection.method = 'mean.var.plot',
                                      mean.cutoff = mean.cutoff,
                                      dispersion.cutoff = dispersion.cutoff)
    seuratObj <- markStepRun(seuratObj, 'FindVariableFeatures', saveFile)
  }

  if (!hasStepRun(seuratObj, 'ScaleData')) {
    seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = c("nCount_RNA", "percent.mito"))
    seuratObj <- markStepRun(seuratObj, 'ScaleData', saveFile)
  }

  if (!hasStepRun(seuratObj, 'RunPCA')) {
    seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(object = seuratObj), verbose = FALSE)
    seuratObj <- markStepRun(seuratObj, 'RunPCA', saveFile)
  }


  if (!hasStepRun(seuratObj, 'ProjectDim')) {
    seuratObj <- ProjectDim(object = seuratObj)
    seuratObj <- markStepRun(seuratObj, 'ProjectDim', saveFile)
  }

  if (!hasStepRun(seuratObj, 'JackStraw')) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100)
    seuratObj <- markStepRun(seuratObj, 'JackStraw', saveFile)
  }

  if (!hasStepRun(seuratObj, 'ScoreJackStraw')) {
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- markStepRun(seuratObj, 'ScoreJackStraw')
  }

  return(seuratObj)
}

PreProcess_SerObjs <- function(SerObj.path = NULL, SerObjRDSKey="SeuratObj.rds",
                               ProjName="10X",
                               save.path = NULL, save.fig.path = NULL,
                               returnList=F,
                               save.fig = T,
                               ENSMB.tag="ENSMM", RhesusConvDavid = F,
                               nUMI.high = 20000, nGene.high = 3000, pMito.high = 0.15,
                               nUMI.low = 0.99, nGene.low = 200, pMito.low = -Inf,
                               RhesusConvDavid.path = "./data/Rhesus/David6.8_ConvertedRhesus_ENSMMUG.txt",
                               fvg.x.low.cutoff = 0.01, fvg.x.high.cutoff = 4.5, fvg.y.cutoff = 1.5,
                               KeepGene.LS =NULL,
                               nDimPCA=15, RemoveCellCycle=F, path2CCfiles="./data/CellCycle", cleanGeneNames=T){



  SerObj.path = "/Volumes/Maggie/Work/OHSU/Bimber/Expts/214/10X/bin/all";
  ProjName="10X";
  save.path = NULL
  returnList = F
  save.fig = T
  save.fig.path = NULL
  ENSMB.tag="ENSMM"
  RhesusConvDavid = T
  nUMI.high = 20000
  nGene.high = 3000
  pMito.high = 0.15
  nUMI.low = 0.99
  nGene.low = 200
  pMito.low = -Inf
  RhesusConvDavid.path = "/Volumes/Maggie/Work/OHSU/Eisa/R/scRNASeq/data/Rhesus/David6.8_ConvertedRhesus_ENSMMUG.txt"
  fvg.x.low.cutoff = 0.01
  fvg.x.high.cutoff = 4.5
  fvg.y.cutoff = 1.5

  KeepGene.LS =SGS.LS
  nDimPCA=15
  RemoveCellCycle=F
  path2CCfiles = "/Volumes/Maggie/Work/OHSU/Eisa/R/scRNASeq/data/CellCycle"
  SerObjRDSKey="SeuratObj.rds"


  require(Seurat)

  if(returnList) TempLS <- list()

  if(!is.null(SerObj.path)){
    if(is.null(save.path)) {
      save.path <- paste(SerObj.path, "/SerProc", sep="")
      if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    }
    if(is.null(save.fig.path)) {
      save.fig.path <- paste(SerObj.path, "/SerProcFigs", sep="")
      if(!dir.exists(save.fig.path)) dir.create(save.fig.path, recursive = T)
    }



    all_RDS  <- list.files(SerObj.path, full.names = T, pattern = ".rds")
    SeurObj_RDS <-  all_RDS[grep(SerObjRDSKey, all_RDS)]

    print("Found files... examples:")
    print(head(SeurObj_RDS))


    for(xN in 1:length(SeurObj_RDS)){

      if(!file.exists(paste(save.path, "/", basename(SeurObj_RDS[xN]), "_proc.rds", sep=""))){

        # xN=1
        if(returnList) TempLS$SeuratObjs <- list()
        print("Reading in...")
        print(SeurObj_RDS[xN])


        SeuratObjs <- readRDS(SeurObj_RDS[xN])

        ### call fx here....
        processSeurat1_SERIII(SeuratObjs,
                              dispersion.cutoff = c(fvg.y.cutoff, Inf),
                              mean.cutoff = c(fvg.x.low.cutoff, fvg.x.high.cutoff))



        TotalPerGeneExpressed      <- rowSums(SeuratObjs@raw.data)
        TotalPerGeneExpressed.perc <- round(TotalPerGeneExpressed/sum(TotalPerGeneExpressed)*100, 5)

        TotalPerCellExpressed <- colSums(SeuratObjs@raw.data)
        TotalPerCellExpressed.perc <- round(TotalPerCellExpressed/sum(TotalPerCellExpressed)*100, 5)

        SeuratObjs <- AddMetaData(object = SeuratObjs,
                                  metadata = TotalPerCellExpressed,
                                  col.name = "SumTotGeneExprPerCell")
        SeuratObjs <- AddMetaData(object = SeuratObjs,
                                  metadata = TotalPerCellExpressed.perc,
                                  col.name = "PercentSumTotGeneExprPerCell")



        if(save.fig) png(filename =  paste(save.fig.path, "/ViolinPlotTrio_preFilt_", basename(SeurObj_RDS[xN]), ".png", sep=""), width = 15, height = 10, units = "in", res=200)
        VlnPlot(object = SeuratObjs,
                features.plot = c("nGene", "nUMI", "percent.mito", "PercentSumTotGeneExprPerCell"),
                nCol = 3, cols.use = col_vector, x.lab.rot=T, size.x.use = 11)
        if(save.fig) dev.off()


        #Genes that dont map to a specific name
        noGeneSYM <- rownames(SeuratObjs@raw.data)[grepl(ENSMB.tag, rownames(SeuratObjs@raw.data))]

        length(noGeneSYM)

        # write.table(noGeneSYM,
        #             "./10X/Rhesus_ENSMMUG.csv",
        #             sep=", ", , row.names = F, quote = F,
        #             col.names = F)


        if(RhesusConvDavid){
          print("Reading in David Data...")

          David6.8ConvTable <- data.frame(read.csv(RhesusConvDavid.path, sep = "\t", header = T))
          rownames(David6.8ConvTable) <- David6.8ConvTable$From
          David6.8ConvTable <- David6.8ConvTable[noGeneSYM, ]
          length(unique(noGeneSYM)); length((noGeneSYM))
          rownames(David6.8ConvTable) <- noGeneSYM

          David6.8ConvTable$Final <- as.character(David6.8ConvTable$To)

          David6.8ConvTable$Final[which(is.na(David6.8ConvTable$To))] <- rownames(David6.8ConvTable)[which(is.na(David6.8ConvTable$To))]

          rownames(SeuratObjs@raw.data)[grepl("ENSMM", rownames(SeuratObjs@raw.data))] <- David6.8ConvTable$Final

          length(rownames(SeuratObjs@raw.data)); length(unique(rownames(SeuratObjs@raw.data)))


          duplicatedGeneNames <- names(which(table(rownames(SeuratObjs@raw.data))>1))

          #change the second duplace name to add a .2
          #perhaps can avg the expr?
          for(geneN in duplicatedGeneNames){
            rownames(SeuratObjs@raw.data)[which(rownames(SeuratObjs@raw.data)==geneN)[2]] <- paste(geneN, ".2", sep="")

          }

          rownames(SeuratObjs@data) <- rownames(SeuratObjs@raw.data)
        }


        print("Filtering Cells...")

        SeuratObjs <- FilterCells(object = SeuratObjs,
                                  subset.names = c("nUMI", "nGene", "percent.mito"),
                                  low.thresholds = c(nUMI.low,   nGene.low,     pMito.low),
                                  high.thresholds = c(nUMI.high, nGene.high,    pMito.high))

        if(save.fig) png(filename =  paste(save.fig.path, "/ViolinPlotTrio_postFilt.png", sep=""), width = 15, height = 10, units = "in", res=200)
        VlnPlot(object = SeuratObjs,
                features.plot = c("nGene", "nUMI", "percent.mito", "PercentSumTotGeneExprPerCell"),
                nCol = 3, cols.use = col_vector, x.lab.rot=T, size.x.use = 11)
        #CleaningLS$Figs$Combo$ViolinPlotTrio_postFilt <- recordPlot()
        if(save.fig) dev.off()

        print("Normalizing ...")

        SeuratObjs <- NormalizeData(object = SeuratObjs,
                                    normalization.method = "LogNormalize",
                                    scale.factor = 10000)



        print("Scaling, centering, and regressing out nUMI and p.mito ...")
        SeuratObjs <- ScaleData(object = SeuratObjs, vars.to.regress = c("nUMI", "percent.mito"))

        print("Finding Variable Genes ...")
        SeuratObjs <- FindVariableGenes(object = SeuratObjs,
                                        mean.function = ExpMean,
                                        dispersion.function = LogVMR,
                                        x.low.cutoff = fvg.x.low.cutoff, #X-axis function is the mean expression level
                                        x.high.cutoff = fvg.x.high.cutoff, #based on plot viz
                                        y.cutoff = fvg.y.cutoff, #Y-axis it is the log(Variance/mean)
                                        num.bin = 20) #y.cutoff = 1 sd away from averge within a bin


        if(!is.null(KeepGene.LS)){
          print("updated Var genes with additional set")
          length(SeuratObjs@var.genes)
          SeuratObjs@var.genes <- unique(c(SeuratObjs@var.genes, as.character(unlist(KeepGene.LS))))
          length(SeuratObjs@var.genes)
        }


        print("Running PCA with Var genes ...")
        SeuratObjs <- RunPCA(object = SeuratObjs,
                             pc.genes = SeuratObjs@var.genes,
                             do.print = F)

        SeuratObjs <- ProjectPCA(object = SeuratObjs, do.print = FALSE)


        if(RemoveCellCycle){
          print("performing cell cycle cleaning ...")

          cc.genes <- readLines(con = paste(path2CCfiles, "/regev_lab_cell_cycle_genes.txt", sep=""))
          g2m.genes <- readLines(con =  paste(path2CCfiles, "/G2M.txt", sep=""))

          # We can segregate this list into markers of G2/M phase and markers of S
          # phase
          s.genes <- cc.genes[1:43]
          g2m.genes <- unique(c(g2m.genes, cc.genes[44:97]))

          s.genes <- s.genes[which(s.genes %in% rownames(SeuratObjs@raw.data))]
          g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(SeuratObjs@raw.data))]



          #OrigPCA <- SeuratObjs@dr$pca
          print("running PCA with cell cycle genes")
          SeuratObjs <- RunPCA(object = SeuratObjs, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
          SeuratObjs <- ProjectPCA(object = SeuratObjs, do.print = FALSE)

          SeuratObjs <- CellCycleScoring(object = SeuratObjs,
                                         s.genes = s.genes,
                                         g2m.genes = g2m.genes,
                                         set.ident = TRUE)

          SeuratObjsCCPCA <- as.data.frame(SeuratObjs@dr$pca@cell.embeddings)
          colnames(SeuratObjsCCPCA) <- paste(colnames(SeuratObjsCCPCA), "CellCycle", sep="_")
          #add the PCA DF to meta.data of SeuratObjs

          if(save.fig) png(filename =  paste(save.fig.path, "/PCAplot_CellCycle_", basename(SeurObj_RDS[xN]), ".png", sep=""), width = 10, height = 10, units = "in", res=200)
          PCAPlot(object = SeuratObjs, dim.1 = 1, dim.2 = 2)
          if(save.fig) dev.off()

          print("regressing out S and G2M score ...")
          SeuratObjs <- ScaleData(object = SeuratObjs,
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  display.progress = T)



          print("Running PCA with Var genes ...")
          SeuratObjs <- RunPCA(object = SeuratObjs,
                               pc.genes = SeuratObjs@var.genes,
                               do.print = F)

          SeuratObjs <- ProjectPCA(object = SeuratObjs, do.print = FALSE)

          SeuratObjs@meta.data <- as.data.frame(cbind(SeuratObjs@meta.data, SeuratObjsCCPCA[rownames(SeuratObjs@meta.data),]))



        }






        if(save.fig) png(filename =  paste(save.fig.path, "/PCElbowPlot_", basename(SeurObj_RDS[xN]), ".png", sep=""), width = 10, height = 10, units = "in", res=200)
        PCElbowPlot(SeuratObjs)
        if(save.fig) dev.off()

        if(save.fig) png(filename =  paste(save.fig.path, "/PCAHeatmap_", basename(SeurObj_RDS[xN]), ".png", sep=""), width = 10, height = 10, units = "in", res=200)
        PCHeatmap(object = SeuratObjs,
                  pc.use = 1:nDimPCA,
                  cells.use = 200,
                  do.balanced = TRUE,
                  label.columns = FALSE,
                  use.full = FALSE)
        if(save.fig) dev.off()

        print("Clustering in PCA Space  ...")
        SeuratObjs <- FindClusters(object = SeuratObjs,
                                   reduction.type = "pca",
                                   dims.use = 1:nDimPCA,
                                   resolution = 0.6,
                                   print.output = 0,
                                   save.SNN = TRUE)

        SeuratObjs <- StashIdent(object = SeuratObjs,
                                 save.name = "Cluster_PCA_0.6")

        if(save.fig) png(filename =  paste(save.fig.path, "/PCAplot_clust_", basename(SeurObj_RDS[xN]), ".png", sep=""), width = 10, height = 10, units = "in", res=200)
        PCAPlot(object = SeuratObjs, dim.1 = 1, dim.2 = 2)
        if(save.fig) dev.off()


        print("Running TSNE on PCA ...")
        SeuratObjs <- RunTSNE(object = SeuratObjs,
                              reduction.type = "pca",
                              dims.use = 1:nDimPCA)

        if(save.fig) png(filename =  paste(save.fig.path, "/TSNEplot_clust_", basename(SeurObj_RDS[xN]), ".png", sep=""), width = 10, height = 10, units = "in", res=200)
        TSNEPlot(object = SeuratObjs, do.label = TRUE)
        if(save.fig) dev.off()


        print("Running UMAP on PCA...")
        SeuratObjs <- RunUMAP(SeuratObjs,
                              cells.use = NULL,
                              dims.use = 1:nDimPCA,
                              reduction.use = "pca",
                              max.dim = 2L,
                              reduction.name = "umap",
                              reduction.key = "UMAP",
                              n_neighbors = 30L,
                              min_dist = 0.3,
                              metric = "correlation",
                              seed.use = 42)


        if(save.fig) png(filename =  paste(save.fig.path, "/UMAPplot_clust_", basename(SeurObj_RDS[xN]), ".png", sep=""), width = 10, height = 10, units = "in", res=200)
        DimPlot(object = SeuratObjs, reduction.use = 'umap')
        if(save.fig) dev.off()




        print("saving ...")
        saveRDS(SeuratObjs,
                paste(save.path, "/", basename(SeurObj_RDS[xN]), "_proc.rds", sep=""))

        if(returnList) TempLS$SeuratObjs[[basename(exp.dirs[xN])]] <- SeuratObjs

      } else {
        if(returnList) {
          #to return a full list of data, need to read what is already processed and saved...
          TempLS$SeuratObjs[[basename(exp.dirs[xN])]] <- readRDS(
            paste(save.path, "/", basename(SeurObj_RDS[xN]), "_proc.rds", sep=""))
        }

      }




    }

    if(returnList) return(TempLS)

  }
}
