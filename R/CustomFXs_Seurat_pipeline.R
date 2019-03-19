
predict_MCE <- function(ProcSERobj.path = NULL, PatternOfProcSERobj="_proc.rds",
                        classification.path = NULL, file.select = NULL,
                        TrainedClassifiers.path = "../PBMC3k/data",
                        save.fig.path = NULL, col_vector=NULL, returnLS = F, GarnettClassify=F,
                        Garnett.path = "./data/Garnett/pbmc_classification.txt", MCEClassify=T,
                        ModuleScoreGeneListClassify=F, ModuleScoreGeneLists=NULL){


  if(is.null(ProcSERobj.path)){
    print("path does not exists")
  }else {
    ClassifiersLS$MCEyhat <- list()
    ClassifiersLS$MCEyhat$CD8T <- list()
    ClassifiersLS$MCEyhat$CD4T <- list()
    ClassifiersLS$MCEyhat$NK <- list()
    ClassifiersLS$MCEyhat$B <- list()
    ClassifiersLS$MCEyhat$Lymph <- list()
    if(GarnettClassify) ClassifiersLS$Garnett <- list()
    if(ModuleScoreGeneListClassify) ClassifiersLS$SeuratGeneScore <- list()

    SERObjects_processed.paths <- list.files(ProcSERobj.path, full.names = T, pattern = PatternOfProcSERobj)

    if(length(SERObjects_processed.paths)==0){
      print("no files found...")
    }else {

      if(is.null(classification.path)) classification.path <- ProcSERobj.path
      if(is.null(save.fig.path)) save.fig.path <- ProcSERobj.path
      if(is.null(col_vector)) col_vector <- colors(distinct = T)

      if(!is.null(file.select)) SERObjects_processed.paths <- SERObjects_processed.paths[grepl(file.select, SERObjects_processed.paths)]

      for(SERObj.path in SERObjects_processed.paths){
        #SERObj.path = SERObjects_processed.paths[1]
        print(basename(SERObj.path))
        tempSER <- readRDS(SERObj.path)

        tempName <- basename(gsub("_", "", gsub("-", "_", gsub("\\.", "", gsub("_SeuratObj.rds_proc.rds", "", SERObj.path)))))



        if(GarnettClassify) {
          print("Starting Garnett Classifier Pipe...")
          if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Garnettyhat.rds",sep=""))){

          ClassifiersLS$Garnett[[tempName]] <- Garnett_Classification_Seurat(tempSER,
                                        marker_file_path = "./data/Garnett/pbmc_classification.txt",
                                        reutrnMonObj=F)
          print("Saving Garnett Results")
          saveRDS(ClassifiersLS$Garnett[[tempName]],
                  paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Garnettyhat.rds",sep=""))


          } else {
            print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Garnettyhat.rds",sep=""))
            print("Already done ... loading for LS")
            ClassifiersLS$Garnett[[tempName]] <- list(GarC=readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Garnettyhat.rds",sep="")))
          }
        }


        if(ModuleScoreGeneListClassify){

            print("Starting Seurat's AddModule Scoring for GeneSets")
            tempLSScores <- list()

            for(GeneList in names(ModuleScoreGeneLists)){
              print(GeneList)

              tempLSScores[[GeneList]] <- AddModuleScoreV2(object=tempSER,
                               genes.list = list(ModuleScoreGeneLists[[GeneList]]),
                               genes.pool = NULL,
                               n.bin = 25,
                          seed.use = 1,
                          ctrl.size = 100,
                          use.k = FALSE,
                          enrich.name = "Cluster",
                          random.seed = 1, returnSerObj=F)
              colnames(tempLSScores[[GeneList]]) <- GeneList

            }

            #Seurate gene score (SGS)
            ClassifiersLS$SeuratGeneScore[[tempName]] <- list(SGS=tempLSScores)

            ##save as.data.frame(tempLSScores) as tsv/csv for use if needed
            }



        #CD8 T cells

        if(MCEClassify){

          if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep=""))){

            #one can directly give the Seurat object to the ClassifyCellsCustom()
            #since looping, its faster to compute the non-sparse log once

            X.SerObj.temp <- log10(Matrix::as.matrix(t(tempSER@data))+1)

            dim(X.SerObj.temp)
            # head(rownames(X.SerObj.temp))


            ClassifiersLS$MCEyhat$CD8T[[tempName]]  <- list(MCE=ClassifyCellsCustom(
              Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_CD8T.rds", sep=""),
              testing.data = X.SerObj.temp, log10T=F))

            ClassifiersLS$MCEyhat$CD8T[[tempName]]$MCE$log10T = T

            #either save results as tsv or csv
            #ClassifiersLS$MCEyhat$CD8T[[tempName]]$yhat.DF
            #or entire object
            saveRDS(ClassifiersLS$MCEyhat$CD8T[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep=""))


          } else {
            print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep=""))
            print("Already done ... loading for LS")
            ClassifiersLS$MCEyhat$CD8T[[tempName]] <- list(MCE=readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD8T_MCEyhat.rds",sep="")))
          }




          #CD4T cells
          tempName <- basename(gsub("_", "", gsub("-", "_", gsub("\\.", "", gsub("_SeuratObj.rds_proc.rds", "", SERObj.path)))))

          if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep=""))){

            if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

            ClassifiersLS$MCEyhat$CD4T[[tempName]]  <- list(MCE=ClassifyCellsCustom(
              Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_CD4T.rds", sep=""),
              testing.data = X.SerObj.temp, log10T=F))
            ClassifiersLS$MCEyhat$CD4T[[tempName]]$MCE$log10T = T

            #either save results as tsv or csv
            #ClassifiersLS$MCEyhat$CD4T[[tempName]]$yhat.DF
            #or entire object
            saveRDS(ClassifiersLS$MCEyhat$CD4T[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep=""))


          } else {
            print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep=""))
            print("Already done ... loading for LS")
            ClassifiersLS$MCEyhat$CD4T[[tempName]] <- list(MCE=readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_CD4T_MCEyhat.rds",sep="")))
          }


          #NK cells
          tempName <- basename(gsub("_", "", gsub("-", "_", gsub("\\.", "", gsub("_SeuratObj.rds_proc.rds", "", SERObj.path)))))

          if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep=""))){

            if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

            ClassifiersLS$MCEyhat$NK[[tempName]]  <- list(MCE=ClassifyCellsCustom(
              Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_NK.rds", sep=""),
              testing.data = X.SerObj.temp, log10T=F))
            ClassifiersLS$MCEyhat$NK[[tempName]]$MCE$log10T = T

            #either save results as tsv or csv
            #ClassifiersLS$MCEyhat$NK[[tempName]]$yhat.DF
            #or entire object
            saveRDS(ClassifiersLS$MCEyhat$NK[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep=""))


          } else {
            print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep=""))
            print("Already done ... loading for LS")
            ClassifiersLS$MCEyhat$NK[[tempName]] <- list(MCE=readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_NK_MCEyhat.rds",sep="")))
          }



          #B cells
          tempName <- basename(gsub("_", "", gsub("-", "_", gsub("\\.", "", gsub("_SeuratObj.rds_proc.rds", "", SERObj.path)))))

          if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep=""))){

            if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

            ClassifiersLS$MCEyhat$B[[tempName]]  <- list(MCE=ClassifyCellsCustom(
              Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_B.rds", sep=""),
              testing.data = X.SerObj.temp, log10T=F))
            ClassifiersLS$MCEyhat$B[[tempName]]$MCE$log10T = T

            #either save results as tsv or csv
            #ClassifiersLS$MCEyhat$B[[tempName]]$yhat.DF
            #or entire object
            saveRDS(ClassifiersLS$MCEyhat$B[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep=""))


          } else {
            print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep=""))
            print("Already done ... loading for LS")
            ClassifiersLS$MCEyhat$B[[tempName]] <- list(MCE=readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_B_MCEyhat.rds",sep="")))
          }



          #Lymph cells
          tempName <- basename(gsub("_", "", gsub("-", "_", gsub("\\.", "", gsub("_SeuratObj.rds_proc.rds", "", SERObj.path)))))

          if(!file.exists(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep=""))){

            if(!exists("X.SerObj.temp")) X.SerObj.temp <- readRDS(SERObj.path)

            ClassifiersLS$MCEyhat$Lymph[[tempName]]  <- list(MCE=ClassifyCellsCustom(
              Classifier.rds.path = paste(TrainedClassifiers.path, "/MCR_LS_Lymph.rds", sep=""),
              testing.data = X.SerObj.temp, log10T=F))

            ClassifiersLS$MCEyhat$Lymph[[tempName]]$MCE$log10T = T

            #either save results as tsv or csv
            #ClassifiersLS$MCEyhat$Lymph[[tempName]]$yhat.DF
            #or entire object
            saveRDS(ClassifiersLS$MCEyhat$Lymph[[tempName]], paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep=""))
          } else {
            print(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep=""))
            print("Already done ... loading for LS")
            ClassifiersLS$MCEyhat$Lymph[[tempName]] <- list(MCE=readRDS(paste(classification.path, "/", basename(gsub(".rds_proc.rds", "", SERObj.path)), "_Lymph_MCEyhat.rds",sep="")))
          }

          print("All classifiers are done with this file ... ")
          print("Saving figures...")

          png(filename = paste(save.fig.path, "/UMAP_PentaClass_", basename(gsub(".rds_proc.rds", "", SERObj.path)), ".png", sep=""), width = 13, height = 9, units = "in", res=200)

          grid.arrange(grobs=list(ggUMAP(tempSER, colFac = cut(round((1-ClassifiersLS$MCEyhat$CD8T[[tempName]]$MCE$yhat.DF$NotProb), 2),breaks=c(-Inf,.5, 0.8, 1)), col_vector = col_vector, ptSize = .9, add2title = "\nPredicted Prob. CD8T"),
                                  ggUMAP(tempSER, colFac = cut(round((1-ClassifiersLS$MCEyhat$CD4T[[tempName]]$MCE$yhat.DF$NotProb), 2),breaks=c(-Inf,.5, 0.8, 1)), col_vector = col_vector, ptSize = .9, add2title = "\nPredicted Prob. CD4T"),
                                  ggUMAP(tempSER, colFac = cut(round((1-ClassifiersLS$MCEyhat$NK[[tempName]]$MCE$yhat.DF$NotProb), 2),breaks=c(-Inf,.5, 0.8, 1)), col_vector = col_vector, ptSize = .9, add2title = "\nPredicted Prob. NK"),
                                  ggUMAP(tempSER, colFac = cut(round((1-ClassifiersLS$MCEyhat$B[[tempName]]$MCE$yhat.DF$NotProb), 2),breaks=c(-Inf,.5, 0.8, 1)), col_vector = col_vector, ptSize = .9, add2title = "\nPredicted Prob. B"),
                                  ggUMAP(tempSER, colFac = cut(round((1-ClassifiersLS$MCEyhat$Lymph[[tempName]]$MCE$yhat.DF$NotProb), 2),breaks=c(-Inf,.5, 0.8, 1)), col_vector = col_vector, ptSize = .9, add2title = "\nPredicted Prob. Lymph")),
                       nrow=2, heights = c(1,1))
          dev.off()


          yhat.Combo <- as.data.frame(cbind(1-ClassifiersLS$MCEyhat$CD8T[[tempName]]$MCE$yhat.DF$NotProb,
                                            1-ClassifiersLS$MCEyhat$CD4T[[tempName]]$MCE$yhat.DF$NotProb,
                                            1-ClassifiersLS$MCEyhat$NK[[tempName]]$MCE$yhat.DF$NotProb,
                                            1-ClassifiersLS$MCEyhat$B[[tempName]]$MCE$yhat.DF$NotProb,
                                            1-ClassifiersLS$MCEyhat$Lymph[[tempName]]$MCE$yhat.DF$NotProb))

          colnames(yhat.Combo) <- c("CD8T", "CD4T", "NK", "B", "Lymph")
          rownames(yhat.Combo) <- colnames(tempSER@data)

          Classificatio.meta.data <- list()

          Classificatio.meta.data$CD8T_MCE <- apply(yhat.Combo, 1, function(xR){
            ifelse(xR["Lymph"] > 0.5 &
                     xR["CD8T"] > 0.5 &
                     xR["CD4T"] < 0.5 &
                     xR["NK"] < 0.5 &
                     xR["B"] < 0.5, 1, 0)
          })

          Classificatio.meta.data$CD4T_MCE <- apply(yhat.Combo, 1, function(xR){
            ifelse(xR["Lymph"] > 0.5 &
                     xR["CD8T"] < 0.5 &
                     xR["CD4T"] > 0.5 &
                     xR["NK"] < 0.5 &
                     xR["B"] < 0.5, 1, 0)
          })


          Classificatio.meta.data$B_MCE <- apply(yhat.Combo, 1, function(xR){
            ifelse(xR["Lymph"] > 0.5 &
                     xR["CD8T"] < 0.5 &
                     xR["CD4T"] < 0.5 &
                     xR["NK"] < 0.5 &
                     xR["B"] > 0.5, 1, 0)
          })


          Classificatio.meta.data$NK_MCE <- apply(yhat.Combo, 1, function(xR){
            ifelse(xR["Lymph"] > 0.5 &
                     xR["CD8T"] < 0.5 &
                     xR["CD4T"] < 0.5 &
                     xR["NK"] > 0.5 &
                     xR["B"] < 0.5, 1, 0)
          })


          Classificatio.meta.data$LymphNotTBNK_MCE <- apply(yhat.Combo, 1, function(xR){
            ifelse(xR["Lymph"] > 0.5 &
                     xR["CD8T"] < 0.5 &
                     xR["CD4T"] < 0.5 &
                     xR["NK"] < 0.5 &
                     xR["B"] < 0.5, 1, 0)
          })

          Classificatio.meta.data$NotLymphTBNK_MCE <- apply(yhat.Combo, 1, function(xR){
            ifelse(xR["Lymph"] < 0.5 &
                     xR["CD8T"] < 0.5 &
                     xR["CD4T"] < 0.5 &
                     xR["NK"] < 0.5 &
                     xR["B"] < 0.5, 1, 0)
          })

          Classificatio.meta.data <- as.data.frame(Classificatio.meta.data)



          png(filename = paste(save.fig.path, "/UMAP_PentaClass_", basename(gsub(".rds_proc.rds", "", SERObj.path)), ".png", sep=""), width = 13, height = 9, units = "in", res=200)

          grid.arrange(grobs=list(ggUMAP(tempSER,
                                         colFac = NULL,
                                         cells.use = rownames(subset(Classificatio.meta.data, NotLymphTBNK_MCE == 1)),
                                         ptSize = .7, ptAlpha = .8,
                                         col_vector = col_vector, add2title="\nNon-Lymphocytes, Non-B-T-NK\nEnsmble Classifier"),
                                  ggUMAP(tempSER,
                                         colFac = NULL,
                                         cells.use = rownames(subset(Classificatio.meta.data, LymphNotTBNK_MCE == 1)),
                                         ptSize = .7, ptAlpha = .8,
                                         col_vector = col_vector, add2title="\nLymphocytes, Non-B-T-NK\nEnsmble Classifier"),
                                  ggUMAP(tempSER,
                                         colFac = NULL,
                                         cells.use = rownames(subset(Classificatio.meta.data, NK_MCE == 1)),
                                         ptSize = .7, ptAlpha = .8,
                                         col_vector = col_vector, add2title="\nNK Lymphocytes, Non-B-T\nEnsmble Classifier"),
                                  ggUMAP(tempSER,
                                         colFac = NULL,
                                         cells.use = rownames(subset(Classificatio.meta.data, B_MCE == 1)),
                                         ptSize = .7, ptAlpha = .8,
                                         col_vector = col_vector, add2title="\nB Lymphocytes, Non-NK-T\nEnsmble Classifier"),
                                  ggUMAP(tempSER,
                                         colFac = NULL,
                                         cells.use = rownames(subset(Classificatio.meta.data, CD4T_MCE == 1)),
                                         ptSize = .7, ptAlpha = .8,
                                         col_vector = col_vector, add2title="\nCD4T Lymphocytes, Non-NK-B-CD8T\nEnsmble Classifier"),
                                  ggUMAP(tempSER,
                                         colFac = NULL,
                                         cells.use = rownames(subset(Classificatio.meta.data, CD8T_MCE == 1)),
                                         ptSize = .7, ptAlpha = .8,
                                         col_vector = col_vector, add2title="\nCD8T Lymphocytes, Non-NK-B-CD4T\nEnsmble Classifier")),
                       nrow=2, heights = c(1,1))
          dev.off()




          ClassifiersLS$classificationLS[[tempName]] <- list(Classificatio.meta.data = Classificatio.meta.data,
                                                             yhat.Combo=yhat.Combo)


        }



      }



    }
    if(returnLS) return(ClassifiersLS)

  }





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



        mito.genes <- grep(pattern = "^MT-", x = rownames(x = SeuratObjs@raw.data), value = TRUE)
        length(mito.genes)

        percent.mito <- Matrix::colSums(SeuratObjs@raw.data[mito.genes, ]) / Matrix::colSums(SeuratObjs@raw.data)


        SeuratObjs <- AddMetaData(object = SeuratObjs,
                                  metadata = percent.mito,
                                  col.name = "percent.mito")

        if(cleanGeneNames){
          rownames(SeuratObjs@data)      <- gsub("-", "", gsub("_", "", rownames(SeuratObjs@data)))
          rownames(SeuratObjs@raw.data)  <- gsub("-", "", gsub("_", "", rownames(SeuratObjs@raw.data)))
        }


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
