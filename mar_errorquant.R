# # script to quantify human error in the marsupial dataset for senstivity analyses
# based on von Cramon-Taubadel et al 2007, Lockwood et al 2002, Fruciano 2016

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

#explore to get a sense of what you expect to find. 'taxa$rep' is a column identifying replicates.
# 'taxa$taxon' is a column identifying different species groups
plot(PCA[[60]]$x[,1:2],bg=c("black","blue")[taxa$rep+1],pch=c(21,22)[taxa$taxon])
anderson(PCA[[47]]$sdev^2) 

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

# # # Addition post-review: evaluating error through ANOVA
# #code development, testing. Note that if you use PCA$x then results identical to those from $scaled
# testgdf<-geomorph.data.frame(coords=PCA[[1]]$scaled,specimen=taxa$specnum)
# testgdf2<-geomorph.data.frame(coords=PCA[[1]]$x,specimen=taxa$specnum)
# errorANOVA<-procD.lm(coords~factor(specimen),data=testgdf,iter=999,RRPP=TRUE) %>% .$aov.table
# repeatability<-((errorANOVA$MS[1]-errorANOVA$MS[2])/6)/(errorANOVA$MS[2]+((errorANOVA$MS[1]-errorANOVA$MS[2])/6))

# test for repeatibility in the dataset as a whole
repvals<-genRepeatvals(PCA,cluster,variable=taxa$specnum,rep=6,pcs=length(groups))
summary_stats<-Rvalsumm(repvals) #organize results
row.names(summary_stats)<-cluster #calculate and label summary statistics for procustes ANOVA scores
write.csv(summary_stats,file="humerr_ProANOVA_summary-stats.csv")

# make figure
tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="humerr_repeatability-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab="repeatability",cex.axis=cex.axis,
          legend.pos='top',legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()

# test for repeatibility PC by PC
identity<-makegroups(PCA,cluster[c(13:16)])
repeatability<-sapply(cluster[c(13:16)],function(x) NULL)
for (cls in 1:length(identity)){
  repeatability[[cls]]<-sapply(identity[[cls]],function(x) NULL) #pairwise correlations of PC1's
  for (i in 1:length(repeatability[[cls]])){
    repeatability[[cls]][[i]]<-find_repeatablePCs(PCA[[identity[[cls]][i]]]$x,variable=taxa$specnum,rep=6)
  }
  repeatability[[cls]]<-unlist(repeatability[[cls]]) %>% 
    matrix(.,nrow=length(repeatability[[cls]]),byrow=TRUE) %>% t
}

repeatability.mat<-unlist(repeatability[[1]]) %>% matrix(.,nrow=nrow(taxa),byrow=FALSE)
for (i in 2:length(repeatability)){
  repeatability.mat<-unlist(repeatability[[i]]) %>% matrix(.,nrow=nrow(taxa),byrow=FALSE) %>%
    rbind(repeatability.mat,.)
}

summary_stats<-lapply(seq_len(nrow(repeatability.mat)), function(i) unlist(repeatability.mat[i,])) %>%
  Rvalsumm 

write.csv(summary_stats,"humerr_ProANOVA_byPC_summary-stats.csv")


# make figure
tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="humerr_repeatPC-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm=seq(1,nrow(taxa)),summary_stats[,4],summary_stats[,5],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Principal Components",ylab="repeatability",cex.axis=cex.axis,
          legend.pos='topright',legend.txt=c("128","256","512","1,024"),legend.title="Pseudolandmarks",
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()

### If you use this code in published materials, please cite: 
# Vitek, N.S., Manz, C.L., Gao, T. Bloch, J.I., Strait, S.G., Boyer, D.M. In Press. Semi-supervised determination of pseudocryptic morphotypes using observer-free characterizations of anatomical alignment and shape. Ecology and Evolution. 
