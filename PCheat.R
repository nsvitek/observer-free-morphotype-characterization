### plot heat map differences along PCs
#input: input data (output from preprocess function), replicate choice (#), palette)
#       whether to alter differences (changes color distribution; alter="square","cube","none")
#       note: if multiple alignments, probably scaled=data
#output: two colored 3D surfaces in ply format, saved.
PCheat<-function(data,scaled,choice=NA,pc=1,palette,alter="none"){
  if (is.na(choice)){
    meanshape<-mshp(scaled$m2d)
    pc.differences<-pcdif(data,meanshape,pcs=pc)
  }
  else{
    meanshape<-mshp(scaled[[choice]]$m2d)
    pc.differences<-pcdif(data[[choice]],meanshape,pcs=pc)
  }
  heatcolor.index<-shpdif(pc.differences[[1]]$min,pc.differences[[1]]$max,palette,alter=alter)
  pieces2plot<-list(pc.differences,heatcolor.index)
  return(pieces2plot)
}