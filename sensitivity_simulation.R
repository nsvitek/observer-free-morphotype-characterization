locateScripts<-"C:/cygwin/home/N.S/scripts/observer-free-morphotype-characterization/"
# location of data
# locateData<-"C:/Users/N.S/Dropbox/Documents/Dissertation/sensitivity_analysis/data"
locateData<-"D:/Dropbox/Documents/Dissertation/sensitivity_analysis/data"

# Load Dependencies, Common Objects ------------------------------------------------------------------
setwd(locateScripts)
source("sensitivity_dependencies.R")
filename<-"morphologika_unscaled_high.txt"
groups<-c(8,1,5,6,7,8,8,8,8,8,1,1,1,1,1,2,3,4) #set group ID's
sim.palette<-brewer.pal(n=8,"Spectral") #set color palette

cluster<-c("0128","0256","0512","1024","2048","4096") #group replicates
# plot settings
pseudolm.lab<-c(128,256,512,1024,2048,4096) #build grouping vectors
col.tab.discrete<-brewer.pal((length(cluster)/length(pseudolm.lab)),"Set2")
ylab.txt<-parse(text="R^2")
cex=1.5
cex.lab=1
cex.axis=1
mtext.line=2
mtext.cex=1
line.lwd=0.2
legend.pos='bottomright'
legend.cex=1.5
pset<-21

taxon.point<-add.alpha(sim.palette,alpha=1)
taxon.bubble<-add.alpha(sim.palette,alpha=0.3)
palette<-colorRampPalette(c("blue","green","yellow","red"))

# Load Data ------------------------------------------------------------------
setwd(locateData)
setwd("simulation")
PCA<-readrep(19,22,c(1:length(groups)),filename=filename) #read in data

# lsdr<-list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
# setwd(lsdr[32])
# objs<-read.morphologika(filename)
# scld<-preprocess(objs)
# shps<-prcomp(scld$m2d,scale.=FALSE)
# anderson(shps$sdev)
# # cbind(dimnames(objs)[[3]],groups) #check to make sure labels are correct

# Analyses for Publication --------------------------------------------------
# # # start analyzing
# PCA[[1]]$x[,1]
# summary(lm(PCA[[1]]$x[,1]~PCA[[4]]$x[,1]))$r.squared
setwd("../outputr")
(r_vals<-genRvals(PCA,cluster))
summary_stats<-Rvalsumm(r_vals)
row.names(summary_stats)<-cluster
write.csv(summary_stats,file="sim_r-vals_summary-stats.csv")
pfish<-Pvalsumm(r_vals,cluster,metric="mean")
write.csv(pfish,"sim_r_pairwise_mean_p-vals.csv",quote=F)
pfish<-Pvalsumm(r_vals,cluster,metric="median")
write.csv(pfish,"sim_r_pairwise_median_p-vals.csv",quote=F)

# pick the two groups with highest and lowest mean R^2, look at distributions
minset<-r_vals[[which(summary_stats[,1]==min(summary_stats[,1]))]]
maxset<-r_vals[[which(summary_stats[,1]==max(summary_stats[,1]))]]
plotRdistr(minset,maxset)

# plot line graph of mean R^2
tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="sim_r-val-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab=ylab.txt,cex.axis=cex.axis,
          mtext.line=mtext.line)
dev.off()

# plot first two PCs with alignment error
choice<-1
point.set<-pset
point.index<-groups
point.color<-taxon.point
bubble.color<-taxon.bubble
cex<-cex
cex.lab<-cex.lab
mtext.line<-3

alignerrPC(PCA,cluster,choice=6,pcs=c(1,2), 
           point.set=rep(21,8),point.index=groups,point.color=sim.palette,
           bubble.color=taxon.bubble,cex=cex,cex.lab=cex.lab,mtext.line=3)

alignerrPC(PCA,cluster,choice=6,pcs=c(1,2),
           point.set=pset,point.index=groups,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=cex.lab,mtext.line=3)

tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="sim_align_1.tif")
par(mar=c(4,5,.5,.5)) #best
alignerrPC(PCA,cluster,choice=5,pcs=c(1,2),
              point.set=rep(21,8),point.index=groups,point.color=taxon.point,
              bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="sim_align_2.tif")
par(mar=c(4,5,.5,.5)) #worst
alignerrPC(PCA,cluster,choice=1,pcs=c(1,2),
              point.set=rep(21,8),point.index=groups,point.color=taxon.point,
              bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="sim_align_3.tif")
par(mar=c(4,5,.5,.5)) #middle
alignerrPC(PCA,cluster,choice=4,pcs=c(1,2),
           point.set=rep(21,8),point.index=groups,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
legend('center',legend=c("shape 1","shape 8"),pch=21,cex=3,
       pt.bg=sim.palette[c(1,8)],pt.cex=6)

# PC Plots ------------------------------------------------------------------
plot(shps$x[,1:2],pch=21,bg=sim.palette[groups],cex=1.5,main=lsdr[32])
text(shps$x[,1],shps$x[,2],dimnames(objs)[[3]],pos=4)
plot3d(PCA[[37]]$x[,1:3],col=sim.palette[groups],size=10)
plot(PCA[[37]]$x[,1:2],pch=21,bg=sim.palette[groups],cex=1.5)

# View PC1 Shapes -------------------------------------------------
diff1<-PCheat(shps,scld,pc=1,palette=palette,alter="none")
open3d()
plot3d(diff1[[1]][[1]]$max,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")
open3d()
plot3d(diff1[[1]][[1]]$min,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")
diff2<-PCheat(shps,scld,pc=2,palette=palette,alter="square")
open3d()
plot3d(diff2[[1]][[1]]$max,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")
open3d()
plot3d(diff2[[1]][[1]]$min,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")

# VarianceDisparity -------------------------------------------------
# are more or less complex (# of patches, size of patches) shapes
# any harder to align properly, i.e., is their more variance among some identical shapes 
# than others?
# From Zelditch et al. 2012 workbook (361-362): In studies of shape, a variance can be calculated by measuring 
# the Procrustes distance of each individual from the mean, which is equivalent to measuring the 
# variance of each coordinate, summed over all the coordinates. Unlike P.D.P.'s function below,
# Workbook divides by N-1, not N

#for 1 & 8, take which(group=i), calculate disparity for that group n terms of PC scores
group1<-which(groups==1) #sphere
group2<-which(groups==8) #water molecule

disparity<-matrix(NA, nrow=2,ncol=length(PCA))
rownames(disparity)<-c("sphere1","water8")
for (i in 1:length(PCA)){
  disparity[1,i]<-individual.disparity(PCA[[i]]$x[group1,])
  disparity[2,i]<-individual.disparity(PCA[[i]]$x[group2,])
}

disparitydiff<-disparity[,] %>% #use only "good" alignments
  apply(.,2, function(x) x[1]/x[2]) %>% #ratio of variance/disparity
  cbind(.,c(rep(128,9),rep(256,9),rep(512,9),rep(1024,9),rep(2048,9),rep(4096,9))) %>% as.data.frame
colnames(disparitydiff)<-c("variance.ratio","id")

# library(ggplot2)
ggplot(data=disparitydiff, aes(x=variance.ratio,fill=factor(id))) +
  geom_histogram(binwidth=0.25) #ratio of disparity depends on # of points

# Note in example below that code from P.D.P and Zelditch et al. calculations produce same results
# individual.disparity(PCA[[37]]$x[group1,])
# apply(PCA[[37]]$m2d[group1,],2,var) %>% sum(.)/(nrow(d2)-1)


# Mantel tests
# mantel_vals<-genMantelvals(PCA,cluster) #takes forever. save results and never overwrite.
# mantelR<-unlist(mantel_vals[[1]]) %>% matrix(.,ncol=36,byrow=TRUE)
# rownames(mantelR)<-cluster
# write.csv(mantelR,"sim_mantelR.csv")
# mantelP<-unlist(mantel_vals[[2]]) %>% matrix(.,ncol=36,byrow=TRUE)
# rownames(mantelP)<-cluster
# write.csv(mantelP,"sim_mantelP.csv")
# summary_stats<-Rvalsumm(mantel_vals[[1]])
mantelR<-read.csv("sim_mantelR.csv",row.names=1,header=TRUE)
summary_stats<-lapply(seq_len(nrow(mantelR)), function(i) unlist(mantelR[i,])) %>%
  Rvalsumm
row.names(summary_stats)<-cluster
# write.csv(summary_stats,file="sim_mantelR_summary-stats.csv")
# pfish<-Pvalsumm(mantel_vals[[1]],cluster,metric="mean")
# write.csv(pfish,"sim_r_pairwise_mean_mantelR.csv",quote=F)
# pfish<-Pvalsumm(mantel_vals[[1]],cluster,metric="median") 
# write.csv(pfish,"sim_r_pairwise_median_mantelR.csv",quote=F)

tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="sim_mantelR-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab="mean observed correlation coefficeint",cex.axis=cex.axis,
          legend.pos='topleft',legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()

#Phenetic trees
RFdists<-getRFdist(PCA,cluster)
unlist(RFdists) %>% matrix(.,ncol=36,byrow=TRUE) %>% 
  write.csv("sim_RFdists_raw.csv")
summary_stats<-Rvalsumm(RFdists)
row.names(summary_stats)<-cluster
write.csv(summary_stats,file="sim_RFdists_summary-stats.csv")
pfish<-Pvalsumm(RFdists,cluster,metric="mean")
write.csv(pfish,"sim_r_pairwise_mean_RFdists.csv",quote=F)
pfish<-Pvalsumm(RFdists,cluster,metric="median")
write.csv(pfish,"sim_r_pairwise_median_RFdists.csv",quote=F)

tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="sim_RFdist-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab="Robinson-Foulds Distance",cex.axis=cex.axis,
          legend.pos='topright',legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()
