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
cl3<-clusterExperiment(assay(cl),clMatNew,transformation=transformation(cl))

test_that("`plotClusters` works with matrix, ClusterExperiment objects", {
    #test matrix version
    x<-plotClusters(clusters=allClusters(cl))
    expect_equal(dim(allClusters(cl)),dim(x$colors))
    expect_equal(dim(allClusters(cl)),dim(x$aligned))
    expect_equal(length(x$clusterColors),ncol(allClusters(cl)))
    xx<-plotClusters(cl2,orderClusters="clusterMany")
    expect_equal(plotClusters(clusters=cl2,orderClusters="pipeline"),xx)
    
    #test CE version
    x<-plotClusters(cl)
    expect_is(x,"ClusterExperiment")
    expect_equal( x,cl)
    plotClusters(cl,resetOrderSamples=TRUE,resetColors=TRUE) 
    x2<-plotClusters(cl,existingColors="all")
    expect_false(isTRUE(all.equal(x1,x2)))
    #test -1
    plotClusters(cl3) 
    x1<-plotClusters(clusters=cl2,orderClusters="pipeline")

  
  wh<-c(3,4,NCOL(allClusters(cl)))
  x2<-plotClusters(cl,orderClusters=wh)
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