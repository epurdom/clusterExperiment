ps<-c(5,10,50)
cl <- clusterMany(simData,dimReduce="PCA",nPCADims=ps,
                     clusterMethod="pam",ks=2:4,findBestK=c(TRUE,FALSE))
cl2<-addClusters(cl,sample(2:5,size=NCOL(simData),replace=TRUE),type="User")
clMatNew<-apply(allClusters(cl),2,function(x){
	wh<-sample(1:nSamples(cl),size=10)
	x[wh]<- -1
	wh<-sample(1:nSamples(cl),size=10)
	x[wh]<- -2
	return(x)
	})
	#make a new object with -1 values
cl3<-clusterCells(assay(cl),clMatNew,transformation=transformation(cl))

test_that("`plotClusters` works with matrix, ClusterCells objects", {
  plotClusters(cl) 
  plotClusters(cl3) #test -1
  x<-plotClusters(cl2,orderClusters="pipeline")
  expect_equal(dim(allClusters(cl)),dim(x$colors))
  expect_equal(dim(allClusters(cl)),dim(x$aligned))
  expect_equal(length(x$groupToColorLegend),ncol(allClusters(cl)))
  xx<-plotClusters(cl2,orderClusters="clusterMany")
  expect_equal(x,xx)
  wh<-c(3,4,NCOL(allClusters(cl)))
  x2<-plotClusters(cl,orderClusters=wh)
  expect_equal(dim(x2$colors),c(NCOL(simData),3))
  expect_equal(dim(x2$aligned),c(NCOL(simData),3))
  expect_equal(length(x2$groupToColorLegend),3)
  x3<-plotClusters(cl,orderClusters=wh)
  expect_equal(x3,x2)
  x4<-plotClusters(cl,orderClusters=wh[c(3,2,1)])
  expect_false(isTRUE(all.equal(x3,x4)))
  x5<-plotClusters(cl,orderClusters=2)
  plotClusters(cl,metaData=sample(2:5,size=NCOL(simData),replace=TRUE))
  plotClusters(cl,metaData=cbind(sample(2:5,size=NCOL(simData),replace=TRUE),sample(2:5,size=NCOL(simData),replace=TRUE)))
})
#' )
#' colnames(cl$clMat) 
#' #make names shorter for plotting
#' colnames(cl$clMat)<-gsub("TRUE","T",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("FALSE","F",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("k=NA,","",colnames(cl$clMat))
#' par(mar=c(2,10,1,1))
#' out<-plotClusters(cl$clMat,axisLine=-2)
#' out$groupToColorLegend[1:2]
#' head(out$color[out$index,1:2])