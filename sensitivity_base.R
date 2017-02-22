# spine for sensitivity analysis code. 
# Includes repeatedly used objects and figure configurations in sensitivity analysis.
# runs mar_repanalysis, mar_errorquant, mar_repfigs, eul_repanalysis, eul_repfigs

# location of scripts
locateScripts<-"C:/cygwin/home/N.S/scripts/observer-free-morphotype-characterization/"
# location of data
# locateData<-"D:/Dropbox/Documents/Dissertation/sensitivity_analysis/data"
locateData<-"C:/Users/N.S/Dropbox/Documents/Dissertation/sensitivity_analysis/data"


# Load Dependencies ------------------------------------------------------------------
setwd(locateScripts)
source("sensitivity_dependencies.R")

# Commonly Used Objects ------------------------------------------------------------------
cluster<-c("005k0128","005k0256","005k0512","005k1024","010k0128","010k0256",
           "010k0512","010k1024","050k0128","050k0256","050k0512","050k1024",
           "100k0128","100k0256","100k0512","100k1024","fulk0128","fulk0256",
           "fulk0512","fulk1024") #group replicates


# plot settings
pseudolm.lab<-c(128,256,512,1024) #build grouping vectors
col.tab.discrete<-brewer.pal((length(cluster)/length(pseudolm.lab)),"Set2")
legend.txt<-c("Downsampled (5,000)","Downsampled (10,000)",
              "Downsampled (50,000)","Downsampled (100,000)","Not Downsampled")
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

# Make Legend ------------------------------------------------------------------
legend_image <- as.raster(matrix(palette(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
# text(x=1.5, y = c(0,1), labels = c("least","most"))
rasterImage(legend_image, 1, 1, 0,0)

# Error ------------------------------------------------------------------
# take a snapshot of objects in list -- used later
freeze<-ls()

#Human error quantification:
setwd(locateScripts)
source("mar_errorquant.R")
rm(list = setdiff(ls(),freeze)) #clean up environment, also removes freeze

# Marsupials ------------------------------------------------------------------
freeze2<-ls() #NEED THIS?
setwd(locateData)
taxa<-read.csv("mardat.csv",header=TRUE,sep=",")
taxa2<-taxa[which(taxa$rep==0),] #remove replicates
setwd("marsup")
PCA<-readrep(10,18,c(1:3)) #read in data

# hist(r_vals$'100k0128')

# # # start analyzing
setwd("../outputr")
r_vals<-genRvals(PCA,cluster)
summary_stats<-Rvalsumm(r_vals)
row.names(summary_stats)<-cluster
write.csv(summary_stats,file="mar_r-vals_summary-stats.csv")
pfish<-Pvalsumm(r_vals,cluster,metric="mean")
write.csv(pfish,"mar_r_pairwise_mean_p-vals.csv",quote=F)
pfish<-Pvalsumm(r_vals,cluster,metric="median")
write.csv(pfish,"mar_r_pairwise_median_p-vals.csv",quote=F)

# pick the two groups with highest and lowest mean R^2, look at distributions
minset<-r_vals[[which(summary_stats[,1]==min(summary_stats[,1]))]]
maxset<-r_vals[[which(summary_stats[,1]==max(summary_stats[,1]))]]
plotRdistr(minset,maxset)

# plot line graph of mean R^2
tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="mar_r-val-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab=ylab.txt,cex.axis=cex.axis,
          legend.pos=legend.pos,legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()

# plot first two PCs with alignment error
alignerrPC(PCA,cluster,choice=13,pcs=c(1,2),
           point.set=pset,point.index=taxa2$taxon,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=cex.lab,mtext.line=3)

#1,5,6,10,13,17 no good, 2,3,4,7,8,9,11,14,16,18,19,20 ok, 15 needed pc3 but no error
summary_stats[c(1,5,6,10,13,17),1]
summary_stats[c(2,3,4,7,8,9,11,14,15,16,18,19,20),1]
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="mar_align_3.tif")
par(mar=c(4,5,.5,.5))
alignerrPC(PCA,cluster,choice=16,pcs=c(1,2),
           point.set=pset,point.index=taxa2$taxon,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="mar_align_2.tif")
par(mar=c(4,5,.5,.5))
alignerrPC(PCA,cluster,choice=5,pcs=c(1,2),
           point.set=pset,point.index=taxa2$taxon,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="mar_align_1.tif")
par(mar=c(4,5,.5,.5))
alignerrPC(PCA,cluster,choice=9,pcs=c(1,2),
           point.set=pset,point.index=taxa2$taxon,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()

# plot PC shape 
cluster[16] #alternate 50k512 and 10k1024
surchoose<-makegroups(PCA,cluster)[[16]][1]
diff1<-PCheat(PCA,PCA,surchoose,pc=1,palette=palette,alter="square")
open3d()
plot3d(diff1[[1]][[1]]$max,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("mar_PC1_max.ply",format="ascii",pointRadius=0.005)
rgl.close()
open3d()
plot3d(diff1[[1]][[1]]$min,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("mar_PC1_min.ply",format="ascii",pointRadius=0.005)
rgl.close()
diff2<-PCheat(PCA,PCA,surchoose,pc=2,palette=palette,alter="square")
open3d()
plot3d(diff2[[1]][[1]]$max,axes=F,col=diff2[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("mar_PC2_max.ply",format="ascii",pointRadius=0.005)
rgl.close()
open3d()
plot3d(diff2[[1]][[1]]$min,axes=F,col=diff2[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("mar_PC2_min.ply",format="ascii",pointRadius=0.005)
rgl.close()



# clean up environment
rm(list = setdiff(ls(),freeze2))

# Erinaceomorphs ------------------------------------------------------------------
freeze3<-ls() #NEED THIS?
setwd(locateData)
taxa<-read.csv("erinacdat.csv",header=TRUE,sep=",")
setwd("erinac")
PCA<-readrep(12,22,c(1:3)) #read in data

# # # # start analyzing, calculate summary statistics
setwd("../outputr")
r_vals<-genRvals(PCA,cluster)
summary_stats<-Rvalsumm(r_vals) 
row.names(summary_stats)<-cluster
write.csv(summary_stats,file="erin_r-vals_summary-stats.csv")
pfish<-Pvalsumm(r_vals,cluster,metric="mean")
write.csv(pfish,"erin_r_pairwise_mean_p-vals.csv",quote=F)
pfish<-Pvalsumm(r_vals,cluster,metric="median")
write.csv(pfish,"erin_r_pairwise_median_p-vals.csv",quote=F)

# pick the two groups with highest and lowest mean R^2
minset<-r_vals[[which(summary_stats[,1]==min(summary_stats[,1]))]]
maxset<-r_vals[[12]]
plotRdistr(minset,maxset)

# plot line graph
tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="eul_r-val-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab=ylab.txt,
          cex.axis=cex.axis,
          legend.pos=legend.pos,legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex)
dev.off()

# plot first two PCs with alignment error
alignerrPC(PCA,cluster,choice=8,pcs=c(1,2),
           point.set=pset,point.index=taxa$posthoc1,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=cex.lab,mtext.line=3)

alignerrPC(PCA,cluster,choice=16,pcs=c(1,3),
           point.set=pset,point.index=taxa$posthoc1,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=cex.lab,mtext.line=3)

plot3d(PCA[[makegroups(PCA,cluster)[[16]][1]]]$x[,1:3],col=taxon.point[taxa$posthoc1],size=6)

cluster[12]
#in  2d:
#1,2,5,6,7,8,9,10,13,14,15,16,17,18,19,20 no good, 3,4,11,12 ok 
#[15,16 no err, but mixed area, which disappears if you add PC3]
#1,5 no good, 2,6 debatable (in 3d)  ok, #3,4 ok incl pc3
range(summary_stats[c(3,4,11,12,15,16),1])
range(summary_stats[c(1,2,5,6,7,8,9,10,13,14,17,18,19,20),1])

tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="eul_align_3.tif")
par(mar=c(4,5,.5,.5))
alignerrPC(PCA,cluster,choice=8,pcs=c(1,2),
           point.set=pset,point.index=taxa$posthoc1,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="eul_align_2.tif")
par(mar=c(4,5,.5,.5))
alignerrPC(PCA,cluster,choice=1,pcs=c(1,2),
           point.set=pset,point.index=taxa$posthoc1,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="eul_align_1.tif")
par(mar=c(4,5,.5,.5))
alignerrPC(PCA,cluster,choice=12,pcs=c(1,2),
           point.set=pset,point.index=taxa$posthoc1,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()


# based on results above, choose which alignment to use for surface visualization
cluster[15] #or 15 or 12
surchoose<-makegroups(PCA,cluster)[[15]][1]
diff1<-PCheat(PCA,PCA,surchoose,pc=1,palette=palette,alter="square")
open3d()
plot3d(diff1[[1]][[1]]$max,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("eul_PC1_max2.ply",format="ascii",pointRadius=0.005)
rgl.close()
open3d()
plot3d(diff1[[1]][[1]]$min,axes=F,col=diff1[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("eul_PC1_min2.ply",format="ascii",pointRadius=0.005)
rgl.close()
diff2<-PCheat(PCA,PCA,surchoose,pc=2,palette=palette,alter="square")
open3d()
plot3d(diff2[[1]][[1]]$max,axes=F,col=diff2[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("eul_PC2_max2.ply",format="ascii",pointRadius=0.005)
rgl.close()
open3d()
plot3d(diff2[[1]][[1]]$min,axes=F,col=diff2[[2]],size=10,xlab="",ylab="",zlab="")
writePLY("eul_PC2_min2.ply",format="ascii",pointRadius=0.005)
rgl.close()

anderson(PCA[[surchoose]]$sdev)
