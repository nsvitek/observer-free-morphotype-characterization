# functions, objects, etc. unique to sensitivity analyses.


### load, format replicate morphologika files 
#input: position of minimum and maximum character to substring for names, 
#       number of principal components to put into uniform orientation
#output: principal components analysis/output for each  replicate in a list
#make sure: you are working in the directory where replicate directories are located
#optional: can change name of morphologika input files from default.
readrep<-function(charmin, charmax, PCs, filename="morphologika_2_unscaled.txt"){
  lsdr<-list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
  objs<-scld<-shps<-cs_all<-substr(lsdr,charmin,charmax)%>% 
    sapply(function(x) NULL) #make null objects to fill
  for (i in 1:(length(lsdr))){ #loop reads in morphometric data in folders
    setwd(lsdr[i])
    objs[[i]]<-read.morphologika(filename)
    scld[[i]]<-preprocess(objs[[i]])
    shps[[i]]<-prcomp(scld[[i]]$m2d,scale.=FALSE)
    print(paste(i,"of",length(lsdr),"generalized procrustes analyses"))
    setwd("..")
  }
  PCA<-flipPC(shps,PCs) #make sure all PCs are in same orientation
  for (i in 1:length(PCA)){PCA[[i]]<-c(PCA[[i]],scld[[i]])}
  return(PCA)
}

### Organize replicates (PCA) by parameter settings (cluster) in a list
#input: vector of grep-able cluster names that are in the names of replicates
#       loaded using readrep function, plus those readrep-ed replicates
#output: vectors of replicate numbers inside, organized in a list by group
makegroups<-function(PCA,cluster){
  rawlist<-lapply(cluster,grep,names(PCA))
  names(rawlist)<-cluster
  return(rawlist)
}

### Make object of all combinations of replicates within a group
#input: vector of grep-able cluster names that are in the names of replicates
#       loaded using readrep function, plus those readrep-ed replicates
#output: matrix (2 rows) inside, organized in a list by group
makecombo<-function(PCA,cluster){
  rawlist<-makegroups(PCA,cluster)
  combinations<-list()
  for (i in 1:(length(rawlist))){
    combinations[[i]]<-combn(rawlist[[i]],2)
  } #make vectors of combinations of replicates
  names(combinations)<-cluster
  return(combinations)
}

### Calculate R^2 for all combinations of replicates within a group (see makegroups)
#input: vector of grep-able cluster names that are in the names of replicates
#       loaded using readrep function, plus those readrep-ed replicates
#output: vectors of R^2 values for each group, organized in a list by group
genRvals<-function(PCA,cluster){
  combinations<-makecombo(PCA,cluster)
  r_vals<-sapply(cluster,function(x) NULL) #pairwise correlations of PC1's
  for (j in 1:(length(combinations))){
    for (i in 1:(ncol(combinations[[j]]))){
      obj<-lm(PCA[[combinations[[j]][1,i]]]$x[,1]~
                PCA[[combinations[[j]][2,i]]]$x[,1])
      obj2<-sqrt(summary(obj)$r.squared)
      r_vals[[j]][i]<-obj2
    }
  }
  return(r_vals)
}

### Make table of summary statistics for the output of genRvals
#input: output of genRvals
#output: table of mean, standard deviation, median, minimum, and maximum for each group
Rvalsumm<-function(r_vals){
  quartile<-matrix(boxplot(r_vals,plot=F)$stats[c(2,4),],ncol=2,byrow=TRUE)
  summary_stats<-matrix(c(
    unlist(sapply(r_vals,mean)),
    unlist(sapply(r_vals,median)),
    unlist(sapply(r_vals,sd)),
    unlist(sapply(r_vals,min)),
    unlist(sapply(r_vals,max)),quartile),nrow=length(r_vals))
  colnames(summary_stats)<-c("mean","median","sd",
                             "min","max","quartile_low","quartile_high")
  return(summary_stats)
}

### Hijacked pairwise.t.test to run Mann-Whitney U test instead of t-test
#input: same as pairwise.t.test
#output: Bonferroni-corrected p-values in half (triangle) of a symmetric table
pairwise.m.w.test<-function (x, g, p.adjust.method = "bonferroni",...) 
{
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  p.adjust.method <- match.arg(p.adjust.method)
  METHOD <- "pairwise Mann-Whitney U tests"
  compare.levels <- function(i, j) {
    y <- x[which(g==levels(g)[i]|g==levels(g)[j])]
    A <- g[which(g==levels(g)[i]|g==levels(g)[j])]
    test<-wilcox.test(y~A)
    test$p.value
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
              p.adjust.method = p.adjust.method)
  class(ans) <- "pairwise.htest"
  ans
}

### Calculate Bonferroni-corrected p-values for pairwise comparisons of R^2 for groups
#input: output of genRvals
# metric can currently take only values "mean" or "median"
#output: Bonferroni-corrected p-values in half (triangle) of a symmetric table
Pvalsumm<-function(r_vals,cluster,metric = "mean"){
  r_labels<-NULL #look at p-values of pairwise t-tests w/ correction
  for (i in 1:length(r_vals)){r_labels<-c(r_labels,rep(cluster[i],length(r_vals[[i]])))}
  r_labels2<-as.factor(r_labels)
  r_vect<-as.vector(unlist(r_vals))
  r_p<-as.data.frame(cbind(r_vect,r_labels2))
  if (metric == "mean") {
    pfish<-pairwise.t.test(r_p[,1],r_p[,2],p.adjust.method="bonferroni")
  } else if (metric == "median") {
    pfish<-pairwise.m.w.test(r_vect,r_labels,p.adjust.method="bonferroni")
  } else {print('ERROR: metric must equal mean or median')}
  colnames(pfish$p.value)<-levels(r_labels2)[1:(length(levels(r_labels2))-1)]
  rownames(pfish$p.value)<-levels(r_labels2)[2:length(levels(r_labels2))]
  pfish2<-round(pfish$p.value,3) 
  return(pfish2)
}

### plot histogram of distributions of R^2 values for two different groups
#input: two vectors of R^2 values
#output: a layered black and white histogram
plotRdistr<-function(minset,maxset){
  histhi<-hist(maxset,breaks=seq((round(min(minset),2)-.01),1,by=0.005))
  histlo<-hist(minset,breaks=histhi$breaks)
  plot(histhi, col=rgb(0,0,0,1/2), xlim=c(round(min(minset),2),1),
       ylim=c(0,max(histhi$counts)),lty=1,
       xlab=parse(text="R^2"),ylab="",main="")  # first histogram
  plot(histlo, col=rgb(1,1,1,1/2), xlim=c(round(min(minset),2),1),lty=2, add=T)
}

### plot line graph of mean or median R^2 values 
#input: the R^2 values to be plot(input), color lookup table for different lines(colortable), 
#       pseudolandmark values(pseudolm), 
#output: a layered black and white histogram
alignLine<-function(input,colortable,pseudolm,quartile.low,quartile.high,
                    pch=21,cex=1,cex.lab=1,cex.axis=1,xlab="",ylab="",
                    mtext.line=2,mtext.cex=1,line.lwd=2,
                    legend.txt=NULL,legend.pos=NULL,legend.title=NULL,legend.cex=NULL,
                    title.txt="",cex.main=1){
  rp<-length(pseudolm) # number of pseudolandmark settings
  rd<-length(input)/rp # number of downsample settings
  plot(rep(pseudolm,rd),input,bg=colortable[rep(1:rd,each=rp)],
       pch=rep(pset,each=rp),cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,
       xlab="",ylab="",
       xlim=c(min(pseudolm),max(pseudolm)),ylim=c(min(quartile.low),max(quartile.high)))
  mtext(xlab,side=1,line=mtext.line,cex=mtext.cex)
  mtext(ylab,side=2,line=mtext.line,cex=mtext.cex)
  for (i in 1:rd){
    lines(pseudolm,quartile.low[(rp*i-(rp-1)):(rp*i)],lwd=2,
          col=colortable[i])
    lines(pseudolm,quartile.high[(rp*i-(rp-1)):(rp*i)],lwd=2,
          col=colortable[i])
    polygon(c(pseudolm,rev(pseudolm)), 
            c(quartile.low[(rp*i-(rp-1)):(rp*i)],
              rev(quartile.high[(rp*i-(rp-1)):(rp*i)])),
            col=adjustcolor(colortable[i],alpha=0.2), border = NA)
  }
  for (i in 1:rd){
    lines(pseudolm,input[(rp*i-(rp-1)):(rp*i)],col=colortable[i],
          lwd=line.lwd,lty=3)	
  }
  points(rep(pseudolm,rd),input,pch=rep(pset,each=rp),bg=colortable[rep(1:rd,each=rp)],cex=cex)
  if (is.null(legend.txt)==FALSE){
    legend(legend.pos,legend.txt,pch=pch,pt.bg=colortable,pt.cex=legend.cex,
           title=legend.title)
    
  }
   title(title.txt,cex.main=cex.main)
}

### plot heat map differences along PCs
#input: input data (output from preprocess function), replicate choice (#), palette)
#       whether to alter differences (changes color distribution; alter="square","cube","none")
#       note: if multiple alignments, probably scaled=data
#output: two colored 3D surfaces in ply format, saved.
PCheat<-function(data,scaled,choice=NA,pc=1,palette,alter="none"){
  if (is.na(choice)){
    meanshape<-mshp(scaled$m2d)
    pc.differences<-pcdif(data,meanshape,pcs=pc)
  }
  else{
    meanshape<-mshp(scaled[[choice]]$m2d)
    pc.differences<-pcdif(data[[choice]],meanshape,pcs=pc)
  }
  heatcolor.index<-shpdif(pc.differences[[1]]$min,pc.differences[[1]]$max,palette,alter=alter)
  pieces2plot<-list(pc.differences,heatcolor.index)
  return(pieces2plot)
}

### calculate range of PC values a specimen could have due to alignment error
#Input: rawlist, a list of which replicates [numbers] are part of a group
#   shapes, a list of prcomp() outputs for each replicate
#   optional: pcs, a vector of which 2 pcs you want to plot
#   optional: metric: which r function you want to use (default = median absolute deviation )
#Output: an object of values needed to draw error bubbles [draw.ellipse()] in a PC plot
errbvals<-function(groupnum,shps,metric=mad,...){
  if(exists("pcs")==FALSE){pcs<-c(1,2)}
  result<-matrix(data=0,nrow=nrow(shps[[groupnum[1]]]$x),ncol=length(pcs))
  temp<-sapply(pcs,function(x) NULL)
  for (i in 1:length(pcs)){
    temp[[i]]<-shps[[groupnum[1]]]$x[,pcs[i]]
    for (j in 2:length(groupnum)){temp[[i]]<-cbind(temp[[i]],shps[[groupnum[j]]]$x[,pcs[i]])}
  } #pull PC1s, then calculate median absolute deviation (for unpredictable, non-normal distributions)
  for (z in 1:length(pcs)){result[,z]<-apply(temp[[z]],1,function(x) metric(x))} 
  return(result)
}

### draw the error ellipse
#Input: x & y vals for ellipse center, the values returned from errbvals()
#   can also set rgb() vals for col and border
#Output: a drawing layer of error ellipses
errbubble<-function(x,y,errcalc,...){
  if(exists("rgbcol")==FALSE){rgbcol<-rgb(.5,.5,.5,.5)}
  if(exists("rgbbor")==FALSE){rgbbor<-rgb(.5,.5,.5,.5)}
  for (i in 1:nrow(errcalc)){
    draw.ellipse(x[i],y[i],a=errcalc[i,1],
                 b=errcalc[i,2],col=rgbcol,border=rgbbor)	
  }
}

### draw the error ellipse
#Input: multiple-PCA dataset (PCA), clustering vector (cluster), choice of which cluster
#       group to represent (choice), which PCs to plot (pcs), pch values to use (point.set)
#       index to match points to right value (point.index), colors for points(point.color)
#       and their error bubbles (bubble.color), other plotting values
#Output: A PC plot that includes error bubbles
alignerrPC<-function(PCA,cluster,choice,pcs=c(1,2),point.set=21,point.index,point.color,
                     bubble.color,cex=1,cex.lab=1,mtext.line=1){
  group<-makegroups(PCA,cluster)[[choice]]
  index<-as.numeric(point.index)
  and<-anderson(PCA[[group[1]]]$sdev)
  bubbles<-errbvals(group,PCA,metric=mad)
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

### Function 6.5 from Claude (2008), used for Mantel tests
#Input: Two square matrices
#Output: Element-wise correlation coefficient between the lower triangles of each matrix

mantrstat<-function (m1, m2)
{ cor(lower.triang(m1),lower.triang(m2))}


### Function 6.6 from Claude (2008)
#Input: A square matrix (m1), number of landmarks x number of dimensions (n)
#       and number of dimensions (coord)
#Output: permuted matrix

permrowscols<-function (m1,n,coord)
{s <- sample(1:(n/coord))
m1[rep(s,coord), rep(s,coord)]}

### Function  6.7 from Claude (2008)
#Input: Two square matrices (m1 and m2),
#       number of dimensions (coord),
#       number of permutations(nperm), and whether to graph (graph)
#Output: observed statistic (r.stat) and a p value (p)

mantel.t2<-function(m1,m2,coord=1,nperm=1000,graph=FALSE,...)
{n<-nrow(m1)
realz<-mantrstat(m1, m2)
nullstats<-replicate(nperm,mantrstat(m1,permrowscols(m2,n,coord)))
pval <- sum(nullstats > realz)/nperm
if (graph) {
  plot(density(nullstats), type = "l", xlim = c(0,1),...)
  abline(v = realz, col="red")   }
list(r.stat = realz, p = pval)}
