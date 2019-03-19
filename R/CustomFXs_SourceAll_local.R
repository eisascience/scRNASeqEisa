

##### this scripts loads all the custom fxs
Path.to.SF = "/Volumes/Maggie/Work/OHSU/Eisa/R/scRNASeq/rbin"
SourceFiles <- list.files(Path.to.SF, pattern = ".R", full.names = T)[!grepl("SourceAll", list.files(Path.to.SF, 
                                                                                                     pattern = ".R"))]

for(SF in SourceFiles){
  print(SF)
  source(SF)
}
