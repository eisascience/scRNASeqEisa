

MakeSerObjs_10XFolders_SERII <- function(counts.path = NULL,
                                   min.cells = 0,
                                   min.genes = 0,
                                   ProjName="10X",
                                   save.path = NULL,
                                   returnList=F, path.exclude="raw", string.exclude=NULL){

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

        SeuratObjs <- CreateSeuratObject(raw.data = Seurat10X,
                                         min.cells = min.cells,   #genes expressed in >= 5 cells
                                         min.genes = min.genes, #Keep all cells with at least 200 detected genes
                                         project = paste(ProjName, FileNames2Save[xN], sep="_"))


        print("Making matrix sparse...")
        SeuratObjs <- MakeSparse(object = SeuratObjs)

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
