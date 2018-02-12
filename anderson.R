#input: a string of standard deviations for PC's, such as object$sdev from a prcomp'ed dataset. 

#output: eigenvalues, % data explained, and significance of each PC
anderson<-function(sdev){
  eigenvalues<-sdev^2
  percent<-round(eigenvalues/sum(eigenvalues),4)*100
  chisquared1<-NULL
  for (i in 1:length(eigenvalues)){
    arity1<--(length(eigenvalues)-1)*(sum(log(eigenvalues[i:(i+1)])))+(length(eigenvalues)-1)*2*(log(sum(eigenvalues[i:(i+1)])/2))
    chisquared1<-c(chisquared1,arity1)
  }
  and<-round(1-pchisq(chisquared1,2),4)
  result<-list(eigenvalues,percent,and)
  names(result)<-c("eigenvalues","percent","pvals")
  return(result)
}

### If you use this code in published materials, please cite: 
# Vitek, N.S., Manz, C.L., Gao, T. Bloch, J.I., Strait, S.G., Boyer, D.M. In Press. Semi-supervised determination of pseudocryptic morphotypes using observer-free characterizations of anatomical alignment and shape. Ecology and Evolution. 