## Functions to reformat morphologika files for downstream analysis.
#depends: geomorph

# library(geomorph)

#function to change morphologika file format to [n, k*m] matrix
rearrange321<-function(morphologika,m){ #based on function by Aaron Wood
  x<-morphologika[,,1]
  for (i in 2:dim(morphologika)[3]){x<-cbind(x,morphologika[,,i])}
  result<-matrix(0,ncol(x)/m,nrow(x)*m)
  count<-sapply(1:m,function(x) NULL)
  for (j in 1:m){count[[j]]<-seq(j,ncol(result)-(m-j),by=m)
  for (i in 1:nrow(result)){result[i,count[[j]]]<-x[,((i-1)*m)+j]}
  }
  return(result)
}

#want to read morphologika [n,m,k] and write all the info/formatting for downstream
#  1. centroid [#]
#  2. vector of centroid sizes [vector]
#  3. mat2 [n, k*m]; n = specimens, m = dimensions, k = correspondence points
preprocess<-function(morphologika){
  m<-dim(morphologika)[2]
  k<-dim(morphologika)[1]
  n<-dim(morphologika)[3]
  scld<-morphologika
  centroid<-apply(morphologika,2,mean)
  cs<-NULL
  for (h in 1:dim(morphologika)[3]){
    cs[h]<-sqrt(sum((t(t(morphologika[,,h])-centroid))^2))
    scld[,,h]<-morphologika[,,h]/cs[h]
  }
  mat2<-rearrange321(scld,m)
  #return a list with all info, labelled $centroid, $cs, and $m2d[2D matrix], plus $k, $m, $n
  result<-NULL
  result[[1]]<-centroid
  result[[2]]<-cs
  result[[3]]<-mat2
  result[[4]]<-k
  result[[5]]<-m
  result[[6]]<-n
  result[[7]]<-scld
  names(result)<-c("centroid","cs","m2d","k","m","n","scaled")
  return(result)
}

#function to input [n, k*m] matrix and get out a morphologika format.
rearrange123<-function(x,m){
  result<-array(dim=c((ncol(x)/m),m,nrow(x)))
  count<-seq(1,ncol(x),by=m)
  for (j in 1:nrow(x)){
    for (i in 1:(ncol(x)/m)){
      result[i,,j]<-unlist(x[j,(count[i]*(m-(m-1))):((count[i]-1)+m)])	
    }
  }
  return(result)
}

### If you use this code in published materials, please cite: 
# Vitek, N.S., Manz, C.L., Gao, T. Bloch, J.I., Strait, S.G., Boyer, D.M. In Press. Semi-supervised determination of pseudocryptic morphotypes using observer-free characterizations of anatomical alignment and shape. Ecology and Evolution. 