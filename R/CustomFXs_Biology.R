CTL_Immune_GeneList <- function(QuickGO.path="./data/QuickGO"){

  SGS.LS <- list()
  SGS.LS$HighlyActivated <- c("TNFRSF9", "NFKBID", "CCL4", "CCL4L2", "IRF8", "IxNG", "CD83", "CD82", "PLEK", "RGCC") #"ENSMMUG00000013779"
  SGS.LS$LessActivated <- c("LTB", "IL7R", "PTPRCAP", "GZMK")
  SGS.LS$Pheno1 <- c("MT1M", "MGEA5")
  SGS.LS$TCellCanonical <- c("CD3G", "CD3D", "CD3E")
  SGS.LS$CD8Canonical <- c("CD8A", "CD8B")
  SGS.LS$BCellCanonical <- c("CD19", "MS4A1", "CD79A")
  SGS.LS$CD4Canonical <- c("CD4", "IL7R")
  SGS.LS$NKCanonical <- c("GZMK", "LYAR", "NKG7", "GZMA", "GNLY", "FGFBP2", "FCGR3A", "CCL4", "GZMH")

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
