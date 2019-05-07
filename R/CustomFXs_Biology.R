CTL_Immune_GeneList <- function(QuickGO.path="./data/QuickGO"){

  #https://www.rndsystems.com/research-area/lymphoid-lineage-markers

  SGS.LS <- list()
  #### IxNG to IFNG

  SGS.LS$HighlyActivated    <- c("TNFRSF9", "NFKBID", "CCL4", "CCL4L2", "IRF8", "IxNG", "CD83", "CD82", "PLEK", "RGCC") #"ENSMMUG00000013779"
  SGS.LS$LessActivated      <- c("LTB", "IL7R", "PTPRCAP", "GZMK")
  SGS.LS$Pheno1             <- c("MT1M", "MGEA5")
  SGS.LS$TCellCanonical     <- c("CD3G", "CD3D", "CD3E")
  SGS.LS$TCellSecondary     <- c("CD2", "CD5", "CD6")
  SGS.LS$TCellTranscription <- c("MAL", "LAT", "TCRIM", "CD28", "ZAP70", "FYN", "GATA3", "LEF1", "TCF7", "RUNX2", "STAT4", "SATB1")
  SGS.LS$CD8Canonical       <- c("CD8A", "CD8B")
  SGS.LS$CD8Subphenos1      <- c("IL2RB", "KLRC1", "KLRG1", "PRF1", "GNLY", "GZMC", "GZMH", "TBX21")
  SGS.LS$BCellCanonical     <- c("CD19", "MS4A1", "CD79A")
  SGS.LS$BCellSecondary     <- c("CD20", "CD21", "CD22", "FCGR2B") #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1479811/
  SGS.LS$BCellCanonicalV3   <- c("CD24", "CD38", "CD72", "CD74", "CD79B", "CD83", "CD86") #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1479811/
  SGS.LS$CD4Canonical       <- c("CD4", "IL7R")
  SGS.LS$CD4Subphenos1      <- c("NOSIP", "RGS10", "FOXP3", "IL2RA", "LTB")
  SGS.LS$CD4Subphenos2      <- c("ANK3", "MXI1",  "CTSB")
  SGS.LS$NKCanonical        <- c("NCAM1", "CD8B", "GZMK", "LYAR", "NKG7", "GZMA", "GNLY", "FGFBP2", "FCGR3A", "CCL4", "GZMH")
  SGS.LS$Erythrocyte        <- c("HBB", "GYPA", "BCAM", "CD36", "EPO", "HPX", "SLC14A1", "LEPR") #https://www.rndsystems.com/research-area/red-blood-cells--rbcs
  SGS.LS$Megakaryocytes     <- c("PF4", "GP9", "ITGA2B", "TMEM40", "LY6G6F","SEPT5", "PTCRA", "TREML1", "CLDN5", "HGD")
  SGS.LS$Myeloid            <- c("ITGAM", "FUT4", "ANPEP", "CD14", "ITGAX","FCGR3A") #https://www.bio-rad-antibodies.com/flow-myeloid-cell-immunophenotyping.html#1
  SGS.LS$Macrophage         <- c("CD14", "ITGAX","FCGR1A", "CD68", "TFRC", "CD86", "CD163", "TLR2", "TLR4") #https://www.bio-rad-antibodies.com/macrophage-m1-m2-tam-tcr-cd169-cd-markers-antibodies.html
  SGS.LS$PlasmaCells        <- c("SDC1")
  SGS.LS$Eosinophils        <- c("CCR3")
  SGS.LS$Basophil           <- c("LMO4", "ENPP3", "IL3RA")
  SGS.LS$Neutrophil         <- c("CEBPE", "S100A8", "S100A9", "FUT4")
  SGS.LS$Monocytes          <- c("CSF1R", "CD14", "FCGR1A", "CD68", "S100A12", "MS4A7", "CKB", "LILRA3")
  SGS.LS$MonocytesFCGR3A    <- c("HES4", "CDKN1C", "FCGR3A", "MS4A7", "CKB", "LILRA3", "IFITM3", "MS4A4A", "LRRC25")
  SGS.LS$MonocytesCD34p     <- c("CD34", "THY1", "ENG", "KIT", "PROM1" )
  SGS.LS$MAIT               <- c("CD8A", "CD8B", "CD4", "KLRB1", "DPP4")
  SGS.LS$Dendritic          <- c("CD74")
  SGS.LS$WBC                <- c("PTPRC")
  SGS.LS$Granulocyte        <- c("CSF2RB", "CSF3R", "IL1R2", "IL1RN", "IL8RB", "IL13RA1", "FPR1", "MME")
  SGS.LS$Lymphoid           <- c("PTPRC", "CD2", "CD5", "RPL21", "RPL23", "RPL27", "RPL31", "RPL35", "RPS14", "RPS21", "RPS24", "RPS3A", "RPS6", "TGT", "TPI1")






  AnnotationFiles <- list.files(QuickGO.path, pattern = "QuickGO", full.names = T)
  GeneLists <- list()

  for(AnnFile in AnnotationFiles){
    # AnnFile = AnnotationFiles[1]
    tempName = gsub(".txt", "", gsub("_","",gsub("annotations-", "", gsub(".tsv","",gsub("-","",basename(AnnFile))))))
    GeneLists$Extra[[tempName]] <-  read.table(
      AnnFile,
      sep="\t", header=TRUE, row.names = NULL, fill = TRUE )
  }

  SGS.LS$QuickGOgenes <- as.character(data.table::rbindlist(lapply(GeneLists$Extra, function(setX){
    subset(setX, TAXON.ID == 9606)[,c("GO.NAME", "SYMBOL")] #10090 = mouse ; 9606 = human
  }))$SYMBOL)

  return(SGS.LS)
}


RhesusGeneDavidConv <- function(ColNames2Conv, RhesusConvDavid.path, ENSMB.tag = "ENSMM", returnFull=F){

  #Seurat 3 cant update the gene names !!
  # see https://github.com/satijalab/seurat/issues/1207

  print("Reading in David Data...")
  noGeneSYM <- ColNames2Conv[grepl(ENSMB.tag, ColNames2Conv)]


  David6.8ConvTable <- data.frame(read.csv(RhesusConvDavid.path, sep = "\t", header = T))
  rownames(David6.8ConvTable) <- David6.8ConvTable$From
  David6.8ConvTable <- David6.8ConvTable[noGeneSYM, ]
  length(unique(noGeneSYM)); length((noGeneSYM))
  rownames(David6.8ConvTable) <- noGeneSYM

  David6.8ConvTable$Final <- as.character(David6.8ConvTable$To)

  David6.8ConvTable$Final[which(is.na(David6.8ConvTable$To))] <- rownames(David6.8ConvTable)[which(is.na(David6.8ConvTable$To))]

  duplicatedGeneNames <- names(which(table(David6.8ConvTable$Final)>1))

  #change the second duplace name to add a .2
  #perhaps can avg the expr?
  for(geneN in duplicatedGeneNames){

    David6.8ConvTable$Final[grep(geneN , David6.8ConvTable$Final)][2] <- paste(geneN, ".2", sep="")

  }


  if(returnFull) return(David6.8ConvTable) else return(David6.8ConvTable$Final)


}

