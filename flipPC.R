##script to determine whether to reverse axes of PCA of replicate alignments

##overall goal: compare PC1 values of each replicate pairwise to 1st replicate. 
#Along with PC1 values of reverse (-) of PC1. 
#if regular is closer, move on. 
#if reversed is closer, reverse the PC, then move on.
#do the same for PC2. 

#num is which PCs you want flipped. If multiple, must be c(1:2)
flipPC=function(shps,num){
  for (q in min(num):max(num)){
    print(paste("Aligning PC",q))
    for (i in 2:length(shps)){
      asis<-sum((shps[[1]]$x[,q]-shps[[i]]$x[,q])^2)
      flip<-sum((shps[[1]]$x[,q]+shps[[i]]$x[,q])^2)
      if (asis<=flip){
        print(paste("** All is well with alignment",i))
      } else {
        shps[[i]]$x[,q]<-(-shps[[i]]$x[,q])
        print(paste("** Alignment",i,"flipped"))
      }
    }
  }
  return(shps)
}

### If you use this code in published materials, please cite: 
# Vitek, N.S., Manz, C.L., Gao, T. Bloch, J.I., Strait, S.G., Boyer, D.M. In Press. Semi-supervised determination of pseudocryptic morphotypes using observer-free characterizations of anatomical alignment and shape. Ecology and Evolution. 
