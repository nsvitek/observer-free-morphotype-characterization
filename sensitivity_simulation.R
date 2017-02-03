getwd()
library(reshape)
# location of scripts
locateScripts<-"C:/cygwin/home/N.S/scripts/observer-free-morphotype-characterization/"
# location of data
locateData<-"D:/Dropbox/Documents/Dissertation/sensitivity_analysis/data"

# Load Dependencies ------------------------------------------------------------------
setwd(locateScripts)
source("sensitivity_dependencies.R")

# Commonly Used Objects ------------------------------------------------------------------
cluster<-c("005k0128","005k0256","005k0512","005k1024","010k0128","010k0256",
           "010k0512","010k1024","050k0128","050k0256","050k0512","050k1024",
           "100k0128","100k0256","100k0512","100k1024","fulk0128","fulk0256",
           "fulk0512","fulk1024") #group replicates

# Hypothesis of problem ------------------------------------------------------------------
setwd(locateData)
# taxa<-read.csv("mardat.csv",header=TRUE,sep=",") #use marsupials dataset
# taxa2<-taxa[-which(taxa$rep==1),] #remove replicates
# setwd("marsup")
# PCA<-readrep(10,18,c(1:3)) #read in data

taxa<-read.csv("erinacdat.csv",header=TRUE,sep=",") #jk, use erinaceomorphs
setwd("erinac")
PCA<-readrep(12,22,c(1:3)) #read in data

#try a sample alignment from some of the worst setting combinations
testarray<-arrayspecs(PCA[[2]]$m2d,p=PCA[[2]]$k,k=PCA[[2]]$m)
dim(testarray)
testarray[1,,1]
testarray[2,,1]
(testarray[2,,1]-testarray[1,,1])^2 %>% sum() %>% sqrt ->testpt #distance between two points
(testarray[,,1]-testarray[,,2])^2 %>% sum %>% sqrt ->testspec #distance between two specimens

# calculate distance between all pairs of specimens
specimen_distances<-NULL
allspecimens<-combn(1:dim(testarray)[3],2) %>% t() #all combinations of specimens, transposed
for (pair in 1:nrow(allspecimens)){
  answer<-(testarray[,,allspecimens[pair,1]]-testarray[,,allspecimens[pair,2]])^2 %>% sum %>% sqrt
  specimen_distances<-c(specimen_distances,answer)
}

pseudolandmark_distances<-NULL
allpseudolandmarks<-combn(1:dim(testarray)[1],2) %>% t() #all combinations of specimens, transposed
for (pair in 1:nrow(allpseudolandmarks)){
  answer<-(testarray[allpseudolandmarks[pair,1],,1]-testarray[allpseudolandmarks[pair,2],,1])^2 %>% sum() %>% sqrt
  pseudolandmark_distances<-c(pseudolandmark_distances,answer)
}
histlm<-hist(pseudolandmark_distances)
histsp<-hist(specimen_distances)

plot(histlm, col=rgb(0,0,1,1/4), xlim=c(0,0.7),lty=2,xlab="distances",
      main="Comparison of inter-landmark and inter-specimen distances")  # first histogram
plot(histsp, col=rgb(1,0,0,1/4), xlim=c(0,0.7),lty=3, add=T)

#try a sample alignment from some of the best setting combinations
testarray2<-arrayspecs(PCA[[127]]$m2d,p=PCA[[127]]$k,k=PCA[[127]]$m)
(testarray2[2,,1]-testarray2[1,,1])^2 %>% sum() %>% sqrt ->testpt2 #distance between two points
(testarray2[,,1]-testarray2[,,2])^2 %>% sum %>% sqrt ->testspec2 #distance between two specimens

specimen_distances2<-NULL
allspecimens2<-combn(1:dim(testarray2)[3],2) %>% t() #all combinations of specimens, transposed
for (pair in 1:nrow(allspecimens2)){
  answer<-(testarray2[,,allspecimens2[pair,1]]-testarray2[,,allspecimens2[pair,2]])^2 %>% sum %>% sqrt
  specimen_distances2<-c(specimen_distances2,answer)
}

pseudolandmark_distances2<-NULL
allpseudolandmarks2<-combn(1:dim(testarray2)[1],2) %>% t() #all combinations of specimens, transposed
length(allpseudolandmarks2)
length(allspecimens2)
for (pair in 1:nrow(allpseudolandmarks2)){
  answer<-(testarray2[allpseudolandmarks2[pair,1],,1]-testarray2[allpseudolandmarks2[pair,2],,1])^2 %>% sum() %>% sqrt
  pseudolandmark_distances2<-c(pseudolandmark_distances2,answer)
}

histlm2<-hist(pseudolandmark_distances2)
histsp2<-hist(specimen_distances2)

plot(histlm2, col=rgb(0,0,1,1/4), xlim=c(0,0.7),ylim=c(0,100),lty=2,xlab="distances",
      main="Comparison of inter-landmark and inter-specimen distances")  # first histogram
plot(histsp2, col=rgb(1,0,0,1/4), xlim=c(0,0.7),lty=3, add=T)

length(which(pseudolandmark_distances2>=min(specimen_distances2)))/length(pseudolandmark_distances2)
length(which(pseudolandmark_distances>=min(specimen_distances)))/length(pseudolandmark_distances)

# length(which(pseudolandmark_distances>=min(specimen_distances)))/length(pseudolandmark_distances) #marsupial min
# [1] 0.2262549
# > length(which(pseudolandmark_distances>=min(specimen_distances)))/length(pseudolandmark_distances)
# [1] 0.2277313 erinaceomorph minimal example
# > length(which(pseudolandmark_distances>=min(specimen_distances)))/length(pseudolandmark_distances)
# [1] 0.2597195 another erinaceomorph minimal example
# > length(which(pseudolandmark_distances>=min(specimen_distances)))/length(pseudolandmark_distances)
# [1] 0.2217028
# > length(which(pseudolandmark_distances2>=min(specimen_distances2)))/length(pseudolandmark_distances2)
# [1] 0.0954394 a 100k512 example for erinaceomorphs
# > length(which(pseudolandmark_distances2>=min(specimen_distances2)))/length(pseudolandmark_distances2)
# [1] 0.008318441 a 100k1024 marsupial example
makegroups(PCA,cluster)
