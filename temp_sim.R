locateScripts<-"C:/cygwin/home/N.S/scripts"
# location of data
locateData<-"C:/Users/N.S/Dropbox/Documents/Dissertation/sensitivity_analysis/data"

# load dependencies:
setwd(locateScripts)
source("sensitivity_dependencies.R")
filename<-"morphologika_unscaled_high.txt"
groups<-c(1,5,6,7,8,8,8,8,8,8,1,1,1,1,1,2,3,4) #set group ID's
sim.palette<-brewer.pal(n=8,"PiYG") #set color palette

setwd(locateData)
setwd("simulation")
lsdr<-list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
setwd(lsdr[2])
objs<-read.morphologika(filename)
scld<-preprocess(objs)
shps<-prcomp(scld$m2d,scale.=FALSE)
dimnames(objs)[[3]]
cbind(dimnames(objs)[[3]],groups)

plot(shps$x[,1:2],pch=21,bg=sim.palette[groups])
dim(objs)
str(objs)

