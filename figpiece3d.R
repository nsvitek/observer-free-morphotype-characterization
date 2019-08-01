#functions related to figuring geometric morphometric point clouds

### make the mean shape of a group of shapes
#input: matrix where rows = specimens, cols= points; # dimensions in dataset (2,3)
#output: n,m,k,meanshape
#default m=3 (3D)
mshp<-function(x,m=3){
  n<-nrow(x)
  k<-ncol(x)/m
  meanshape<-t(matrix(colMeans(x),m,k))
  result<-list(n,m,k,meanshape)
  names(result)<-c("n","m","k","meanshape")
  return(result)
}

### calculate the shapes at the maxima and minima of principal component axes
#input: object output from prcomp() function, mshp output, desired PCs (ex: 1, 2:4)
#output: k x m shapes for the maximum and minimum of each desired PC
pcdif<-function(prx, mid, pcs=1){
  result<-NULL
  name<-NULL
  j<-1
  for (i in pcs){
    dle<-NULL
    pcmax<-max(prx$rotation[,i])
    pcmin<-min(prx$rotation[,i])
    maxtrans<-t(matrix(pcmax*prx$rotation[,i],mid$m,mid$k))
    mintrans<-t(matrix(pcmin*prx$rotation[,i],mid$m,mid$k))
    maxshp<-mid$meanshape+maxtrans
    minshp<-mid$meanshape+mintrans
    dle$max<-maxshp
    dle$min<-minshp
    result[[j]]<-dle
    name[j]<-paste("pc",i,sep="")
    j<-j+1
  }
  names(result)<-name
  return(result)
}


### calculate the difference between two shapes to express as colors
#input: two shapes in a k x m matrix (row = # pts, col = # dimensions)
#output: color settings for plotting differences between two shapes
#things you could theoretically change: # of colors to cut by,
#includes built-in functions "square" and "cube" change the distribution of distances
shpdif<-function(shape1,shape2,palette,alter="none",outlier=FALSE,...){
  if(exists("ncut")==FALSE){ncut<-256}
  square<-function(x){x^2}
  cube<-function(x){x^3}
  change1<-(shape1-shape2)
  change2<-apply(change1,1,function(x) sqrt(sum(x^2)))
  if (outlier==TRUE){
    change2[which(change2==max(change2))]<-median(change2) #remove outlier
    change2[which(change2==max(change2))]<-median(change2) #remove outlier
    change2[which(change2==max(change2))]<-median(change2) #remove outlier
  }
  if(alter=="none"){col<-change2} else {col<-unlist(lapply(change2,alter))}
  colbydiff<-rescale(col,c(1,ncut)) %>% round
  raintable<-palette(ncut)
  result<-raintable[colbydiff]
  return(result)
}

## Add an alpha value to a colour, from magesblog post 40 april 2013
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

### If you use this code in published materials, please cite: 
# Vitek, N.S., Manz, C.L., Gao, T. Bloch, J.I., Strait, S.G., Boyer, D.M. In Press. Semi-supervised determination of pseudocryptic morphotypes using observer-free characterizations of anatomical alignment and shape. Ecology and Evolution. 