# Mantel Tests ------------------------------------------------------------------
# library(ade4)

#Claude 2008 points out that "Procrustean" data needs a special function
#because rows and columns of landmark data are not independent and interchangeable
#(they are parts of coordinates). He has modified functions accordingly
#resulting in function mantel.t2(), copied to sensitivity_utils.R
#following Claude's (2008) example (p. 265), 
#calculate variance-covariance matrix from 2d array of specimen coordinates
str(PCA[[1]])
dim(PCA[[1]]$m2d)
anderson(PCA[[180]]$sdev)

var(PCA[[180]]$m2d)[1:3,1:3]
mantel.t2(var(PCA[[180]]$m2d),var(PCA[[179]]$m2d),coord=3,nperm=99,graph=TRUE)
test1<-mantel.t2(var(PCA[[1]]$m2d),var(PCA[[4]]$m2d),coord=3,nperm=99,graph=TRUE)
str(test1)

j<-1
i<-3

genMantelvals<-function(PCA,cluster){
  combinations<-makecombo(PCA,cluster)
  mantelR<-sapply(cluster,function(x) NULL) #pairwise correlations of PC1's
  mantelP<-sapply(cluster,function(x) NULL) #pairwise correlations of PC1's
  for (j in 1:1){
    for (i in 1:(ncol(combinations[[j]]))){
      obj<-mantel.t2(var(PCA[[combinations[[j]][1,i]]]$m2d),var(PCA[[combinations[[j]][2,i]]]$m2d),coord=3,nperm=999,graph=FALSE)
      mantelR[[j]][i]<-obj$r.stat
      mantelP[[j]][i]<-obj$p
    }
  }
  return(list(mantelR,mantelP))
}

mantel_vals<-genMantelvals(PCA,cluster)

# Phenograms ------------------------------------------------------------------
library(phangorn)

getRFdist<-function(PCA,cluster){
  combinations<-makecombo(PCA,cluster)
  distRF<-sapply(cluster,function(x) NULL) #pairwise correlations of PC1's
  for (j in 1:1){
    for (i in 1:(ncol(combinations[[j]]))){
      obj<-mantel.t2(var(PCA[[combinations[[j]][1,i]]]$m2d),var(PCA[[combinations[[j]][2,i]]]$m2d),coord=3,nperm=999,graph=FALSE)
      mantelR[[j]][i]<-obj$r.stat
      mantelP[[j]][i]<-obj$p
    }
  }
  return(list(mantelR,mantelP))
  
}
# calculate euclidean distance matrix, 2 ways (coordinates vs. PCs that explain 95% of variation)
for (i in 1:length(anderson(PCA[[180]]$sdev)$percent)){
  enough<-sum(anderson(PCA[[180]]$sdev)$percent[1:i])
  if(enough>=95){
    enough<-i
    break
  }
}
D1<-dist(PCA[[160]]$x[,1:enough], "euclidean") #better?
# perform the clustering;  
cluster<-hclust(D1, method="average") #other option is ward.D, ward.D2 or "average" for UPGMA
# compute the distances along the branches
#in tests with marsupial rep 180, "average" had highest correlation
D2<-cophenetic(cluster)
# compute the correlation and write it to the screen
cor(D1,D2) #<0.85 indicates significant distortion
# draw the dendrogram
cluster$labels<-taxa2$specnum
# plot(cluster,names=taxa2$specnum)
#phangorn version, calculate Robinson-Foulds distance,
#measure of difference between trees. 
#(Steel and Penny 1993; Robinson and Foulds 1981; Wright and Hillis 2014)
distRF[[j]][i]<-RF.dist(as.phylo(cluster),as.phylo(cluster2),check.labels=FALSE)
