locateScripts<-"C:/cygwin/home/N.S/scripts/observer-free-morphotype-characterization/"
# location of data
# locateData<-"C:/Users/N.S/Dropbox/Documents/Dissertation/sensitivity_analysis/data"
locateData<-"D:/Dropbox/Documents/Dissertation/sensitivity_analysis/data"


# load dependencies:
setwd(locateScripts)
source("sensitivity_dependencies.R")
filename<-"morphologika_unscaled_high.txt"
# groups<-c(1,5,6,7,8,8,8,8,8,8,1,1,1,1,1,2,3,4) #set group ID's
groups<-c(rep(1,6),rep(2,6),rep(3,6),rep(4,6),rep(5,6),rep(6,6),rep(7,6),rep(8,6))
sim.palette<-brewer.pal(n=8,"PiYG") #set color palette

setwd(locateData)
setwd("simulation")
lsdr<-list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
setwd(lsdr[5])
objs<-read.morphologika(filename)
scld<-preprocess(objs)
shps<-prcomp(scld$m2d,scale.=FALSE)
cbind(dimnames(objs)[[3]],groups)

plot(shps$x[,1:2],pch=21,bg=sim.palette[groups],cex=1.5)
dim(objs)
str(objs)
anderson(shps$sdev)
str(scld)
plot3d(shps$x[,1:3],col=sim.palette[groups],size=10)
sim.palette
