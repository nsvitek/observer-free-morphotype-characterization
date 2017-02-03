#Mantel Tests
library(ade4)

#covariances function as distance matrices
str(PCA[[1]])
dim(PCA[[1]]$m2d)
anderson(PCA[[1]]$sdev)

PCA[[1]]$x
PCA[[2]]$x

mantel.rtest(, , nrepet = 99)

#Phenograms

#sample code from box turtle, in this case, fos_set is essentially m2d data
D1<-dist(fos_set, "euclidean")
D1<-dist(PCAfos$x[,1:7], "euclidean")
# perform the clustering;  
cluster<-hclust(D1, method="ward.D") #other option is ward.D, ward.D2 or "average" for UPGMA
# compute the distances along the branches
D2<-cophenetic(cluster)
# compute the correlation and write it to the screen
cor(D1,D2) #<0.85 indicates significant distortion
# draw the dendrogram
cluster$labels<-fos_metadata$id
cluster$labels<-fos_metadata$site
plot(cluster,names=fos_metadata$id)