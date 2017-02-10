locateScripts<-"C:/cygwin/home/N.S/scripts/observer-free-morphotype-characterization/"
# location of data
# locateData<-"C:/Users/N.S/Dropbox/Documents/Dissertation/sensitivity_analysis/data"
locateData<-"D:/Dropbox/Documents/Dissertation/sensitivity_analysis/data"

# Load Dependencies, Common Objects ------------------------------------------------------------------
setwd(locateScripts)
source("sensitivity_dependencies.R")
filename<-"morphologika_unscaled_high.txt"
groups<-c(1,5,6,7,8,8,8,8,8,8,1,1,1,1,1,2,3,4) #set group ID's
# groups<-c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(7,5),rep(8,5))
# groups<-c(rep(1,6),rep(2,6),rep(3,6),rep(4,6),rep(5,6),rep(6,6),rep(7,6),rep(8,6))

sim.palette<-brewer.pal(n=8,"PiYG") #set color palette

# Load Data ------------------------------------------------------------------
setwd(locateData)
setwd("simulation2")
lsdr<-list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
setwd(lsdr[5])
objs<-read.morphologika(filename)
scld<-preprocess(objs)
shps<-prcomp(scld$m2d,scale.=FALSE)
# anderson(shps$sdev)
# cbind(dimnames(objs)[[3]],groups) #check to make sure labels are correct

# PC Plots ------------------------------------------------------------------
plot(shps$x[,1:2],pch=21,bg=sim.palette[groups],cex=1.5,main=lsdr[5])
text(shps$x[,1],shps$x[,2],dimnames(objs)[[3]],pos=4)
plot3d(shps$x[,1:3],col=sim.palette[groups],size=10)

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
# any harder to align properly, i.e., is their more variance among some identical shapes than others?
# From Zelditch et al. 2012 workbook (361-362): In studies of shape, a variance can be calculated by measuring 
# the Procrustes distance of each individual from the mean, which is equivalent to measuring the 
# variance of each coordinate, summed over all the coordinates. Unlike P.D.P.'s function below,
# Workbook divides by N-1, not N

############################################################################
#
#   Creates a function called individual.disparity() that calculates the 
#   morphological disparity among several specimens 
#   from continuous data, including PC scores.  Disparity is calculated as 
#   the mean squared distance among the specimens.  There is no standardization
#   because individual specimens do not have variances like groups do.
#
#   The format is: individual.disparity( data )
#
#   where data are the continuous data with individual specimens on
#   rows and variables in columns. Written by P. David Polly, 2008
#
############################################################################


individual.disparity <- function(d) {
  dists <-( dist(d))^2
  return(mean(dists))
}

#for 1:8, take which(group=i), calculate disparity for that group
