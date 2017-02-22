# spine for sensitivity analysis code, Gonzalez et al 2016 dataset. 


# location of scripts
locateScripts<-"C:/Users/N.S/Documents/Dissertation/sensitivity_analysis/code"
# location of data
locateData<-"C:/Users/N.S/Documents/Dissertation/sensitivity_analysis/data"

# load dependencies:
setwd(locateScripts)
source("sensitivity_dependencies.R")

# # # # # commmonly used objects
cluster<-c("ts_0399","ts_0600","ts_1000","ts_2000","ts_3000") #group replicates

# plot settings
pseudolm.lab<-c(399,600,1000,2000,3000) #build grouping vectors
col.tab.discrete<-brewer.pal((length(cluster)/length(pseudolm.lab)),"Set2")
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
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab=ylab.txt,cex.axis=cex.axis,
          legend.pos=legend.pos,legend.txt=legend.txt,legend.title=legend.title,
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

rangetemp<-function(x){
  r<-range(x)
  singlenumber<-r[2]-r[1]
  return(singlenumber)
}
alignerrPCmod<-function(PCA,cluster,choice,pcs=c(1,2),point.set=21,point.index,point.color,
                     bubble.color,cex=1,cex.lab=1,mtext.line=1){
  group<-makegroups(PCA,cluster)[[choice]]
  index<-as.numeric(point.index)
  and<-anderson(PCA[[group[1]]]$sdev)
  bubbles<-errbvals(group,PCA,metric=rangetemp)
  pca<-pcb<-NULL #find median value for each point.
  for (i in 1:length(group)){
    pca<-cbind(pca,PCA[[group[i]]]$x[,pcs[1]])
    pcb<-cbind(pcb,PCA[[group[i]]]$x[,pcs[2]])
  }
  pca.median<-apply(pca,1,median)
  pcb.median<-apply(pcb,1,median)
  toplot<-cbind(pca.median,pcb.median)
  xlim<-c(min(toplot[,pcs[1]]-bubbles[,1]),max(toplot[,pcs[1]]+bubbles[,1]))
  ylim<-c(min(toplot[,pcs[2]]-bubbles[,2]),max(toplot[,pcs[2]]+bubbles[,2]))
  plot(toplot[,pcs[1]],toplot[,pcs[2]],bg=point.color[index],pch=point.set[index],
       cex=cex,xlim=xlim,ylim=ylim,xlab="",ylab="")
  mtext(paste("PC ",pcs[2]," (",round(and$percent[pcs[2]],3),"%)",sep=""),side=2,
        line=mtext.line,cex=cex.lab)
  mtext(paste("PC ",pcs[1]," (",round(and$percent[pcs[1]],3),"%)",sep=""),side=1,
        line=mtext.line,cex=cex.lab)
  for (i in 1:nrow(bubbles)){
    draw.ellipse(toplot[,pcs[1]][i],toplot[,pcs[2]][i],a=bubbles[i,1],b=bubbles[i,2],
                 col=bubble.color[point.index[i]],border=bubble.color[point.index[i]])	
  }
  points(toplot[,pcs[1]],toplot[,pcs[2]],bg=point.color[point.index],
         pch=point.set[point.index],cex=cex)
}

alignerrPCmod(PCA,cluster,choice=1,pcs=c(1,2),
           point.set=pset,point.index=taxa$Experimental.Group,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=cex.lab,mtext.line=3)

alignerrPCmod(PCA,cluster,choice=4,pcs=c(1,2),
              point.set=pset,point.index=taxa$Experimental.Group,point.color=taxon.point,
              bubble.color=taxon.bubble,cex=cex,cex.lab=cex.lab,mtext.line=3)

# plot3d(PCA[[makegroups(PCA,cluster)[[2]][1]]]$x[,1:3],col=taxon.point[taxa$posthoc1],size=6)
plot(PCA[[makegroups(PCA,cluster)[[1]][3]]]$x[,1:2],bg=taxon.point[taxa$Experimental.Group],pch=21)
cluster[12]
#in  2d:
#1,2,5,6,7,8,9,10,13,14,15,16,17,18,19,20 no good, 3,4,11,12 ok 
#[15,16 no err, but mixed area, which disappears if you add PC3]
#1,5 no good, 2,6 debatable (in 3d)  ok, #3,4 ok incl pc3
range(summary_stats[c(3,4,11,12,15,16),1])
range(summary_stats[c(1,2,5,6,7,8,9,10,13,14,17,18,19,20),1])

tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="gon_align_1.tif")
par(mar=c(4,5,.5,.5))
alignerrPCmod(PCA,cluster,choice=1,pcs=c(1,2),
           point.set=pset,point.index=taxa$Experimental.Group,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="gon_align_2.tif")
par(mar=c(4,5,.5,.5))
alignerrPCmod(PCA,cluster,choice=5,pcs=c(1,2),
           point.set=pset,point.index=taxa$Experimental.Group,point.color=taxon.point,
           bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
dev.off()
# tiff(width=7,height=7,units="cm",res=300,pointsize=8,filename="gon_align_1.tif")
# par(mar=c(4,5,.5,.5))
# alignerrPC(PCA,cluster,choice=12,pcs=c(1,2),
#            point.set=pset,point.index=taxa$posthoc1,point.color=taxon.point,
#            bubble.color=taxon.bubble,cex=cex,cex.lab=2,mtext.line=2.5)
# dev.off()


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
