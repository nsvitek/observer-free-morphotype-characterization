##script to determine whether to reverse axes of PCA of replicate alignments

##overall goal: compare PC1 values of each replicate pairwise to 1st replicate. 
#Along with PC1 values of reverse (-) of PC1. 
#if regular is closer, move on. 
#if reversed is closer, reverse the PC, then move on.
#do the same for PC2. 

# library(geomorph)
# setwd("C:/cygwin64/home/N.S/scripts")
# source("rearrange3ddat.R")
# 
# ##Get replicates to work with
# setwd("C:/Users/N.S/Documents/Dissertation/Rwork/mardat_NO_PLY")
# getwd()
# lsdr<-list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
# lsdr<-lsdr[81:89]
# ssr<-substr(lsdr,10,22)#get names, just fulk_128 reps
# objs<-sapply(ssr,function(x) NULL)
# scld<-sapply(ssr,function(x) NULL)
# mat2<-sapply(ssr,function(x) NULL)
# shps<-sapply(ssr,function(x) NULL)
# cs_all<-sapply(ssr,function(x) NULL)
# 
# 
# for (i in 1:(length(lsdr))){ #loop reads in morphometric data in folders
#   setwd(lsdr[i])
#   objs[[i]]<-read.morphologika("morphologika_2_unscaled.txt")
#   shape<-objs[[i]][,,-c(1,11,27)]
#   scld[[i]]<-shape
#   centroid<-apply(shape,2,mean)
#   cs<-NULL
#   for (h in 1:dim(shape)[3]){
#     cs[h]<-sqrt(sum((t(t(shape[,,h])-centroid))^2))
#     scld[[i]][,,h]<-shape[,,h]/cs[h]
#   }
#   cs_all[[i]]<-cs
#   mat1<-scld[[i]][,,1]
#   for (j in 2:dim(scld[[i]])[3]){mat1<-cbind(mat1,scld[[i]][,,j])}
#   mat2[[i]]<-rearrange1(mat1)
#   shps[[i]]<-prcomp(mat2[[i]],scale.=FALSE)
#   print(paste(i,"of",length(lsdr),"generalized procrustes analyses"))
#   setwd("..")
# }

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

