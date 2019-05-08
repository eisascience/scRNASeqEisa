

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

AddModuleScore_SERII <- function(
  object,
  genes.list = NULL,
  genes.pool = NULL,
  n.bin = 25,
  seed.use = 1,
  ctrl.size = 100,
  use.k = FALSE,
  enrich.name = "Cluster",
  random.seed = 1
) {
  set.seed(seed = random.seed)
  genes.old <- genes.list
  if (use.k) {
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(
      X = genes.list,
      FUN = function(x) {
        return(intersect(x = x, y = rownames(x = object@data)))
      }
    )
    cluster.length <- length(x = genes.list)
  }
  if (!all(LengthCheck(values = genes.list))) {
    warning(paste(
      'Could not find enough genes in the object from the following gene lists:',
      paste(names(x = which(x = ! LengthCheck(values = genes.list)))),
      'Attempting to match case...'
    ))
    genes.list <- lapply(
      X = genes.old,
      FUN = CaseMatch, match = rownames(x = object@data)
    )
  }
  if (!all(LengthCheck(values = genes.list))) {
    stop(paste(
      'The following gene lists do not have enough genes present in the object:',
      paste(names(x = which(x = ! LengthCheck(values = genes.list)))),
      'exiting...'
    ))
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(x = object@data)
  }
  data.avg <- Matrix::rowMeans(x = object@data[genes.pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
  ))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
          size = ctrl.size,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object@data)
  )
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = object@data[genes.use, ])
  }
  genes.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object@data)
  )
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    data.use <- object@data[genes.use, , drop = FALSE]
    genes.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
  rownames(x = genes.scores.use) <- colnames(x = object@data)
  object <- AddMetaData(
    object = object,
    metadata = genes.scores.use,
    col.name = colnames(x = genes.scores.use)
  )
  gc(verbose = FALSE)
  return(object)
}
