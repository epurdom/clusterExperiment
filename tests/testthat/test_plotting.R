ps<-c(5,10,50)
cl <- compareChoices(simData,dimReduce="PCA",nPCADims=ps,
   clusterMethod="pam",ks=2:4,findBestK=c(TRUE,FALSE))
test_that("`plotTracking` works with matrix, ClusterCells objects", {
  plotTracking(cl)          
})
#' )
#' colnames(cl$clMat) 
#' #make names shorter for plotting
#' colnames(cl$clMat)<-gsub("TRUE","T",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("FALSE","F",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("k=NA,","",colnames(cl$clMat))
#' par(mar=c(2,10,1,1))
#' out<-plotTracking(cl$clMat,axisLine=-2)
#' out$groupToColorLegend[1:2]
#' head(out$color[out$index,1:2])