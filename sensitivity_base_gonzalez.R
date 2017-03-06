# spine for sensitivity analysis code, Gonzalez et al 2016 dataset. 


# location of scripts
locateScripts<-"C:/cygwin/home/N.S/scripts/observer-free-morphotype-characterization/"
# location of data
# locateData<-"D:/Dropbox/Documents/Dissertation/sensitivity_analysis/data"
locateData<-"C:/Users/N.S/Dropbox/Documents/Dissertation/sensitivity_analysis/data"

# load dependencies:
setwd(locateScripts)
source("sensitivity_dependencies.R")

# # # # # commmonly used objects
cluster<-c("ts_0399","ts_0600","ts_1000","ts_2000","ts_3000") #group replicates

# plot settings
pseudolm.lab<-c(399,600,1000,2000,3000) #build grouping vectors
col.tab.discrete<-brewer.pal(5,"Set2")
# legend.txt<-c("Downsampled (2-4k)","Downsampled (9.3k)")
ylab.txt<-parse(text="R^2")
legend.title<-"Surface Downsampling"
cex=1.5
cex.lab=1
cex.axis=1
mtext.line=2
mtext.cex=1
line.lwd=0.2
legend.pos='bottomright'
legend.cex=1.5
pset<-c(21,22,23,24,25)
taxon.point<-c(rgb(.7,.3,.3,1),rgb(.3,.7,1,1))
taxon.bubble<-c(rgb(.7,.3,.3,.3),rgb(.3,.7,.7,.3))
palette<-colorRampPalette(c("blue","green","yellow","red"))

# # # # # Gonzalez alignment senstivity analysis. 
setwd(locateData)
taxa<-read.csv("gonzalezdat.csv",header=TRUE,sep=",")
#remove replicates
setwd("gonzalez_output")
PCA<-readrep(3,18,c(1:nrow(taxa)),filename="morphologika_unscaled_high.txt") #read in data
# hist(r_vals$'100k0128')

# # # start analyzing
# PCA[[1]]$x[,1]
# summary(lm(PCA[[1]]$x[,1]~PCA[[4]]$x[,1]))$r.squared
setwd("../outputr")
(r_vals<-genRvals(PCA,cluster))
summary_stats<-Rvalsumm(r_vals)
row.names(summary_stats)<-cluster
write.csv(summary_stats,file="gon_r-vals_summary-stats.csv")
pfish<-Pvalsumm(r_vals,cluster,metric="mean")
write.csv(pfish,"gon_r_pairwise_mean_p-vals.csv",quote=F)
pfish<-Pvalsumm(r_vals,cluster,metric="median")
write.csv(pfish,"gon_r_pairwise_median_p-vals.csv",quote=F)

# pick the two groups with highest and lowest mean R^2, look at distributions
minset<-r_vals[[which(summary_stats[,1]==min(summary_stats[,1]))]]
maxset<-r_vals[[which(summary_stats[,1]==max(summary_stats[,1]))]]
plotRdistr(minset,maxset)

# plot line graph of mean R^2
tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="gon_r-val-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete[5],pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset[5],cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab=ylab.txt,cex.axis=cex.axis,
          legend.pos=legend.pos,legend.txt=n,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()

# plot first two PCs with alignment error
choice<-1
point.set<-pset
point.index<-taxa$Experimental.Group
point.color<-taxon.point
bubble.color<-taxon.bubble
cex<-cex
cex.lab<-cex.lab
mtext.line<-3
tst<-c(1,2,4,6,7)
range(tst)

alignerrPC(PCA,cluster,choice=5,pcs=c(1,2),
           point.set=pset,point.index=taxa$Group,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)


# plot3d(PCA[[makegroups(PCA,cluster)[[2]][1]]]$x[,1:3],col=taxon.point[taxa$posthoc1],size=6)
plot(PCA[[makegroups(PCA,cluster)[[1]][3]]]$x[,1:2],bg=taxon.point[taxa$Experimental.Group],pch=21)
cluster[12]
#in  2d:

tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="gon_align_3.tif")
par(mar=c(4,5,.5,.5)) #in between
alignerrPC(PCA,cluster,choice=3,pcs=c(1,2),
           point.set=pset,point.index=taxa$Group,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="gon_align_2.tif")
par(mar=c(4,5,.5,.5)) #best
alignerrPC(PCA,cluster,choice=5,pcs=c(1,2),
           point.set=pset,point.index=taxa$Group,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="gon_align_1.tif")
par(mar=c(4,5,.5,.5)) #worst
alignerrPC(PCA,cluster,choice=1,pcs=c(1,2),
           point.set=pset,point.index=taxa$Group,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()


# based on results above, choose which alignment to use for surface visualization
cluster[5] #or 15 or 12
surchoose<-makegroups(PCA,cluster)[[5]][1]
control<-mshp(PCA[[surchoose]]$m2d[which(taxa$Experimental.Group=="Control"),])
nogh<-mshp(PCA[[surchoose]]$m2d[which(taxa$Experimental.Group=="No GH"),])
colorguide<-shpdif(control$meanshape,nogh$meanshape,palette)

open3d()
plot3d(control$meanshape,axes=F,col=colorguide,size=10,xlab="",ylab="",zlab="")
writePLY("gon_control.ply",format="ascii",pointRadius=0.005)
rgl.close()
open3d()
plot3d(nogh$meanshape,axes=F,col=colorguide,size=10,xlab="",ylab="",zlab="")
writePLY("gon_noGH.ply",format="ascii",pointRadius=0.005)
rgl.close()


# Mantel tests
# mantel_vals<-genMantelvals(PCA,cluster) #takes forever. save results and never overwrite.
# mantelR<-unlist(mantel_vals[[1]]) %>% matrix(.,ncol=36,byrow=TRUE)
# rownames(mantelR)<-cluster
# write.csv(mantelR,"gon_mantelR.csv")
# mantelP<-unlist(mantel_vals[[2]]) %>% matrix(.,ncol=36,byrow=TRUE)
# rownames(mantelP)<-cluster
# write.csv(mantelP,"gon_mantelP.csv")
# summary_stats<-Rvalsumm(mantel_vals[[1]])
mantelR<-read.csv("gon_mantelR.csv",row.names=1,header=TRUE)
summary_stats<-lapply(seq_len(nrow(mantelR)), function(i) unlist(mantelR[i,])) %>%
  Rvalsumm
row.names(summary_stats)<-cluster
# write.csv(summary_stats,file="gon_mantelR_summary-stats.csv")
# pfish<-Pvalsumm(mantel_vals[[1]],cluster,metric="mean")
# write.csv(pfish,"gon_r_pairwise_mean_mantelR.csv",quote=F)
# pfish<-Pvalsumm(mantel_vals[[1]],cluster,metric="median") 
# write.csv(pfish,"gon_r_pairwise_median_mantelR.csv",quote=F)

tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="gon_mantelR-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab="mean observed correlation coefficeint",cex.axis=cex.axis,
          legend.pos='topleft',legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()

#Phenetic trees
RFdists<-getRFdist(PCA,cluster,tips=taxa$Ind)
unlist(RFdists) %>% matrix(.,ncol=36,byrow=TRUE) %>% 
  write.csv("gon_RFdists_raw.csv")
summary_stats<-Rvalsumm(RFdists)
row.names(summary_stats)<-cluster
write.csv(summary_stats,file="gon_RFdists_summary-stats.csv")
pfish<-Pvalsumm(RFdists,cluster,metric="mean")
write.csv(pfish,"gon_r_pairwise_mean_RFdists.csv",quote=F)
pfish<-Pvalsumm(RFdists,cluster,metric="median")
write.csv(pfish,"gon_r_pairwise_median_RFdists.csv",quote=F)

tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="gon_RFdist-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete[5],pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset[5],cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab="Robinson-Foulds Distance",cex.axis=cex.axis,
          legend.pos='n',legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()