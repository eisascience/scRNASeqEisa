
ColorTheme <- function(){
  require(RColorBrewer)
  
  scaleyellowred <- colorRampPalette(c("dodgerblue", "lightyellow", "red"), space = "rgb")(30)
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- col_vector[-4]
  
  return(list(col_vector=col_vector, scaleyellowred=scaleyellowred))
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

MDSmyDF <- function(dfx, labelsDF, factorV, title = "MDS Plot", col_vector){

  
  if(length(factorV) == nrow(dfx)){
    mds <- cmdscale(as.matrix(dist(dfx)))
    colnames(mds) <- c("MDS1", "MDS2")
    
    mds <- cbind(mds, labelsDF) #samples as rows add a cols of labels
    
    
    p1 <- ggplot(mds, aes(x = MDS1, y = MDS2)) +
      theme_bw() +
      geom_hline(yintercept = 0, color = "gray70") +
      geom_vline(xintercept = 0, color = "gray70") +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
      coord_cartesian(xlim = c(min(mds[,1])-5, max(mds[,1])+5)) + 
      scale_color_manual(values = col_vector, name="Samples")
    
    # the graphic with ggrepel
    p1 + geom_text_repel(aes(y = MDS2 + 0.25), label = factorV) +  
      ggtitle(paste("MDS of:",title ,sep=" "))+ 
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    print("rotate t() your dataX")
  }
  
  
}

PCAmyDF <- function (dfx, labels, factorV, title = "PCA Plot", scale, center, col_vector, namePointBL = F) {
  if(class(labels) == "function") {
    print("no labels, using factor as names")
    labels = as.character(factorV)
  }
  if(length(as.character(factorV)) != length(labels)) {
    print("labels length != factor length, using factor as names")
    labels = as.character(factorV)
  }
  
  dfx.pca <- prcomp(t(dfx), scale.=scale, center = center)
  
  MinXY <- min(c(round(min(dfx.pca$rotation[,1:2]) - abs(min(dfx.pca$rotation[,1:2])*0.5),1), -1) )
  MaxXY <- max(c(round(max(dfx.pca$rotation[,1:2]) + abs(max(dfx.pca$rotation[,1:2])*0.5),1),  1) )
  
  
  if(namePointBL){
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) + 
      scale_color_manual(values = col_vector, name="Samples") + 
      geom_text_repel(aes(y = PC2, label = labels))  +  
      ggtitle(paste("PCA of:",title ,sep=" ")) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  } else {
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) + 
      scale_color_manual(values = col_vector, name="Samples")  +  
      ggtitle(paste("PCA of:",title ,sep=" "))+ 
      theme(plot.title = element_text(hjust = 0.5)) + 
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  }
  
  
  
  
}

transposedt <- function(dt, varlabel="myVar") {
  require(data.table)
  dtrows = names(dt)
  dtcols = as.list(c(dt[,1]))
  dtt = transpose(dt)
  dtt[, eval(varlabel) := dtrows]
  setnames(dtt, old = names(dtt), new = c(dtcols[[1]], eval(varlabel)))
  dtt = dtt[-1,]
  setcolorder(dtt, c(eval(varlabel), names(dtt)[1:(ncol(dtt) - 1)]))
  return(dtt)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}





LogAdd <- function(x) {
  mpi <- max(x)
  return(mpi + log(x = sum(exp(x = x - mpi))))
}