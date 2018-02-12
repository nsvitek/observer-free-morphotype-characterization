# Mantel Tests ------------------------------------------------------------------
# library(ade4)

#Claude 2008 points out that Procrustean data needs a special function
#because rows and columns of landmark data are not independent and interchangeable
#(they are parts of coordinates). He has modified functions accordingly
#resulting in function mantel.t2(), copied to sensitivity_utils.R
#following Claude's (2008) example (p. 265), 
#calculate variance-covariance matrix from 2d array of specimen coordinates

genMantelvals<-function(PCA,cluster){
  combinations<-makecombo(PCA,cluster)
  mantelR<-sapply(cluster,function(x) NULL) #pairwise Mantel test R's w/in replicates
  mantelP<-sapply(cluster,function(x) NULL) #significance of R's above
  for (j in 1:(length(combinations))){
    for (i in 1:(ncol(combinations[[j]]))){
      obj<-mantel.t2(var(PCA[[combinations[[j]][1,i]]]$m2d),var(PCA[[combinations[[j]][2,i]]]$m2d),coord=3,nperm=999,graph=FALSE)
      mantelR[[j]][i]<-obj$r.stat
      mantelP[[j]][i]<-obj$p
      print(paste(j,i,sep="."))
    }
  }
  return(list(mantelR,mantelP))
}

mantel_vals<-genMantelvals(PCA,cluster) #takes multiple days to run
mantelR<-unlist(mantel_vals[[1]]) %>% matrix(.,ncol=36,byrow=TRUE)
# mantelR<-unlist(mantel_vals[[1]]) %>% matrix(.,ncol=3,byrow=TRUE)
rownames(mantelR)<-cluster
write.csv(mantelR,"sim_mantelR.csv")
mantelP<-unlist(mantel_vals[[2]]) %>% matrix(.,ncol=36,byrow=TRUE)
# mantelP<-unlist(mantel_vals[[1]]) %>% matrix(.,ncol=3,byrow=TRUE)
rownames(mantelP)<-cluster
write.csv(mantelP,"sim_mantelP.csv")

summary_stats<-Rvalsumm(mantel_vals[[1]])
# mantelR<-read.csv("mar_mantelR.csv",row.names=1,header=TRUE) %>% []
# summary_stats<-lapply(seq_len(nrow(mantelR)), function(i) unlist(mantelR[i,])) %>% 
#   Rvalsumm 
row.names(summary_stats)<-cluster
# write.csv(summary_stats,file="mar_mantelR_summary-stats.csv")
# pfish<-Pvalsumm(mantel_vals[[1]],cluster,metric="mean")
# write.csv(pfish,"mar_r_pairwise_mean_mantelR.csv",quote=F)
# pfish<-Pvalsumm(mantel_vals[[1]],cluster,metric="median")
# write.csv(pfish,"mar_r_pairwise_median_mantelR.csv",quote=F)

tiff(width=7,height=7,units="cm",res=800,pointsize=8,filename="eul_mantelR-mean_line.tif")
par(mar=c(3,3.3,.5,.5))
alignLine(summary_stats[,1],col.tab.discrete,pseudolm.lab,summary_stats[,6],summary_stats[,7],
          pch=pset,cex=cex,cex.lab=cex.lab,xlab="Pseudolandmarks",ylab="mean observed correlation coefficeint",cex.axis=cex.axis,
          legend.pos='topleft',legend.txt=legend.txt,legend.title=legend.title,
          legend.cex=legend.cex,mtext.line=mtext.line)
dev.off()

# Dendrograms ------------------------------------------------------------------
library(phangorn)

getRFdist<-function(PCA,cluster){
  combinations<-makecombo(PCA,cluster)
  distRF<-sapply(cluster,function(x) NULL) #pairwise RF distances of phenograms
  for (j in 1:(length(combinations))){
    for (i in 1:(ncol(combinations[[j]]))){
      for (k in 1:length(anderson(PCA[[combinations[[j]][i]]]$sdev)$percent)){
        enough<-sum(anderson(PCA[[combinations[[j]][i]]]$sdev)$percent[1:k])
        if(enough>=95){
          enough<-k
          break
        }
      }
      D1<-dist(PCA[[combinations[[j]][1,i]]]$x[,1:enough], "euclidean") #better?
      # perform the clustering;  
      cluster1<-hclust(D1, method="average") #other option is ward.D, ward.D2 or "average" for UPGMA
      # compute the distances along the branches
      # D2<-cophenetic(cluster) #for correlation testing, used in development
      # compute the correlation and write it to the screen
      # cor(D1,D2) #<0.85 indicates significant distortion
      # in tests with marsupial rep 180, "average" had highest correlation
      # draw the dendrogram
      cluster2<-dist(PCA[[combinations[[j]][2,i]]]$x[,1:enough], "euclidean") %>%
        hclust(., method="average")
      # cluster2$labels<-cluster1$labels<-taxa$specnum
      # cluster2$labels<-cluster1$labels<-taxa2$specnum
      cluster2$labels<-cluster1$labels<-as.character(groups)
      # plot(cluster,names=taxa2$specnum)
      #phangorn version, calculate Robinson-Foulds distance,measure of difference between trees. 
      #(Steel and Penny 1993; Robinson and Foulds 1981; Wright and Hillis 2014)
      distRF[[j]][i]<-RF.dist(as.phylo(cluster1),as.phylo(cluster2),check.labels=FALSE)
      print(paste(j,i,sep="."))
      }
  }
  return(distRF)
}

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