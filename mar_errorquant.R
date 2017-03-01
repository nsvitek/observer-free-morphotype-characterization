# # script to quantify human error in the marsupial dataset for senstivity analyses
# based on von Cramon-Taubadel et al 2007, Lockwood et al 2002

# sourced by sensitivity_base-objs.R, where things like file locations are stored.

# dependencies: magrittr, vegan, geomorph.
# library(magrittr,vegan)
# library(geomorph)
# setwd(locateScripts)
# source("rearrange3ddat.R")
# source("flipPC.R")
# source("anderson.R")

# load metadata
setwd(locateData) #STARTING FOLDER
taxa<-read.csv("mardat.csv",header=TRUE,sep=",")
# # # # # # # # # # #taxa<-taxa[-c(9,43,40),]
# str(taxa)
repnames<-c("CAB11-0301","CAB13-0887","CAB14-0695") #replicated specimens

# load, format replicate morphologika files 
setwd("marsup-humerr")
PCA<-readrep(10,18,c(1:3))

#start analyzing
setwd("../outputr")

#explore to get a sense of what you expect to find
plot(PCA[[60]]$x[,1:2],bg=c("black","blue")[taxa$rep+1],pch=c(21,22)[taxa$taxon])
anderson(PCA[[47]]$sdev^2) #first four PCs significant

#make a vector of which  specimens are replicated specimens, which are not
repgroups<-rep(1,nrow(taxa)) #base vector of replicate groups, to be modified
for (i in 1:length(repnames)){repgroups[grep(repnames[i],taxa$filename)]<-i+1}

origdist<-errdist<-sapply(PCA,function(x) NULL) #make null objects to fill
for (i in 1:length(PCA)){
  origdist[[i]]<-as.vector(dist(PCA[[i]]$x[which(taxa$rep==0),],method="euclidean"))
  for (j in 1:length(repnames)){
    errdist[[i]]<-c(errdist[[i]],as.vector(dist(PCA[[i]]$x[which(repgroups==(3)),1],method="euclidean")))
  } #calculate pairwise euclidean distances in PC spae between non-replicate specimens
} #then do the same for each replicate specimen individually
#look at distribution of distances
# hist(origdist[[40]]) #spot check: no crazy skewed/bimodal distributions
# hist(unlist(errdist[[40]])) #spot check: no crazy skewed/bimodal distributions
# hist(unlist(sapply(errdist,mean))) #anything crazy here?
# hist(unlist(sapply(origdist,mean)))

# report summary statistics for each  replicate:
# mean, standard deviation, range, ratio of means replicate:non-replicate
summary_stats<-matrix(c(
  unlist(sapply(origdist,mean)),
  unlist(sapply(origdist,sd)),
  unlist(sapply(origdist,min)),
  unlist(sapply(origdist,max)),
  unlist(sapply(errdist,mean)),
  unlist(sapply(errdist,sd)),
  unlist(sapply(errdist,min)),
  unlist(sapply(errdist,max))
),nrow=length(PCA))
summary_stats<-cbind(summary_stats,summary_stats[,5]/summary_stats[,1])
colnames(summary_stats)<-c("original.mean","original.sd","original.min",
                           "original.max","replicate.mean","replicate.sd",
                           "replicate.min","replicate.max","ratio.orig:rep")
head(summary_stats)
write.csv(summary_stats,file="humerr_allreps_summary-stats.csv")

#Also summarize by parameter combination, take mean of each triplicate
rawlist<-lapply(cluster,grep,names(PCA))
summary_groups<-matrix(0,nrow=length(cluster),ncol=ncol(summary_stats))
for (i in 1:length(cluster)){summary_groups[i,]<-apply(summary_stats[rawlist[[i]],],2,mean)}
row.names(summary_groups)<-cluster
colnames(summary_groups)<-colnames(summary_stats)

#write group-wise summary statistics to csv file
write.csv(summary_groups,file="humerr_groups_summary-stats.csv")

#Addition post-review: evaluating error through ANOVA
library(geomorph)

# #code development, testing. Note that if you use PCA$x then results identical to those from $scaled
# testgdf<-geomorph.data.frame(coords=PCA[[1]]$scaled,specimen=taxa$specnum)
# testgdf2<-geomorph.data.frame(coords=PCA[[1]]$x,specimen=taxa$specnum)
# errorANOVA<-procD.lm(coords~factor(specimen),data=testgdf,iter=999,RRPP=TRUE) %>% .$aov.table
# repeatibility<-((errorANOVA$MS[1]-errorANOVA$MS[2])/6)/(errorANOVA$MS[2]+((errorANOVA$MS[1]-errorANOVA$MS[2])/6))

genRepeatvals<-function(PCA,cluster,variable,rep,pcs=length(variable)){
  identity<-makegroups(PCA,cluster)
  repeatibility<-sapply(cluster,function(x) NULL) #pairwise correlations of PC1's
  for (j in 1:(length(identity))){
    for (i in 1:(length(identity[[j]]))){
      testgdf<-geomorph.data.frame(coords=PCA[[identity[[j]][i]]]$x[,pcs],specimen=variable)
      errorANOVA<-procD.lm(coords~factor(specimen),data=testgdf,iter=999,RRPP=TRUE) %>% .$aov.table
      repeatibility[[j]][i]<-((errorANOVA$MS[1]-errorANOVA$MS[2])/rep)/(errorANOVA$MS[2]+((errorANOVA$MS[1]-errorANOVA$MS[2])/rep))
      #6 used in repeatibility because 6 replicates of each repeated specimen
      print(paste(j,i,sep="."))
    }
  }
  return(repeatibility)
}

repeat_vals<-genRepeatvals(PCA,cluster)
repeat_vals2<-genRepeatvals(PCA,cluster)
repeat_vals3<-genRepeatvals(PCA,cluster,variable=taxa$specnum,pcs=1,rep=6)
repeat_vals4<-genRepeatvals(PCA,cluster,variable=taxa$specnum,pcs=2,rep=6)
repeat_vals5<-genRepeatvals(PCA,cluster,variable=taxa$specnum,pcs=3,rep=6)

boxplot(repeat_vals3)
boxplot(repeat_vals4)
boxplot(repeat_vals5)

repeat_val1<-genRepeatvals(PCA,cluster,variable=groups,rep=6,pcs=1) 
repeat_val2<-genRepeatvals(PCA,cluster,variable=groups,rep=6,pcs=2) 
repeat_val3<-genRepeatvals(PCA,cluster,variable=groups,rep=6,pcs=3) 
repeat_val4<-genRepeatvals(PCA,cluster,variable=groups,rep=6,pcs=4) 

boxplot(repeat_val1)
boxplot(repeat_val2)
boxplot(repeat_val3)
boxplot(repeat_val4)

summary_stats<-Rvalsumm(repeat_val1)
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab="repeatibility",cex.axis=cex.axis,
          mtext.line=mtext.line)

#So what kind of summary plot would be helpful for determining how many PCs to use if only one alignment run
# find.repeatablePCs(PCA[[40]]$x,variable=groups,rep=6) for simulated data
step1<-find.repeatablePCs(PCA[[46]]$x,variable=taxa$specnum,rep=6)
step2<-find.repeatablePCs(PCA[[1]]$x,variable=taxa$specnum,rep=6)

makegroups(PCA,cluster)
