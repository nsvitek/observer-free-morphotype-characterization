### generate Procustes ANOVA-based repeatability values for all principal components 
### resulting from a single auto3dgm alignment
#Input: principal components scores from a PCA, for example, PCA$x (PCscores), 
#       how many replicates of a shape are in the dataset (rep), 
#       the vector that labels the replicates (rep)
#Output: vectors of repeatability values for each PC and a scree plot

# library(geomorph) #relies on functions in the geomorph package

find_repeatablePCs<-function(PCscores,variable,rep){
  repeatability<-rep(NA,ncol(PCscores))
  for(i in 1:ncol(PCscores)){
    testgdf<-geomorph.data.frame(coords=PCscores[,i],specimen=factor(variable))
    errorANOVA<-procD.lm(coords~specimen,data=testgdf,iter=999,RRPP=TRUE) %>% .$aov.table
    repeatability[i]<-((errorANOVA$MS[1]-errorANOVA$MS[2])/rep)/(errorANOVA$MS[2]+((errorANOVA$MS[1]-errorANOVA$MS[2])/rep))
  }
  plot(repeatability,xlab="principal components")
  lines(repeatability)
  abline(h=0.95,col="red",lty=2)
  abline(h=0.90,col="blue",lty=3)
  return(repeatability)
}
