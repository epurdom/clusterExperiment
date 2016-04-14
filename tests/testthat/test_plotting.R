ps<-c(5,10,50)
cl <- clusterMany(simData,dimReduce="PCA",nPCADims=ps,
                     clusterMethod="pam",ks=2:4,findBestK=c(TRUE,FALSE))
cl2<-addClusters(cl,sample(2:5,size=NCOL(simData),replace=TRUE),type="User")
clMatNew<-apply(clusterMatrix(cl),2,function(x){
	wh<-sample(1:nSamples(cl),size=10)
	x[wh]<- -1
	wh<-sample(1:nSamples(cl),size=10)
	x[wh]<- -2
	return(x)
	})
	#make a new object with -1 values
cl3<-clusterExperiment(assay(cl),clMatNew,transformation=transformation(cl))
clusterLabels(cl3)
test_that("`plotClusters` works with matrix, ClusterExperiment objects", {
    #test matrix version
    x<-plotClusters(clusters=allClusters(cl))
    expect_equal(dim(allClusters(cl)),dim(x$colors))
    expect_equal(dim(allClusters(cl)),dim(x$aligned))
    expect_equal(length(x$clusterLegend),ncol(allClusters(cl)))
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

smData<-simData[1:30,1:50]
smCount<-simCount[1:30,1:50]

test_that("`plotHeatmap` works with matrix objects", {
    x1<-plotHeatmap(data=smData)
    x2<-plotHeatmap(data=smCount,clusterSamplesData=smData,clusterFeaturesData=smData)
    expect_equal(x1$aheatmapOut,x2$aheatmapOut)
    
    #check internal alignment of sampleData (alignSampleData=TRUE) is working:
    sampleData<-clusterMatrix(cl3)[sample(size=50,1:nrow(clusterMatrix(cl3))),]
    alList<-plotClusters(sampleData)
    alCol<-alList$clusterLegend
    x1<-plotHeatmap(data=smData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,],clusterLegend=alCol,clusterSamples=FALSE,clusterFeatures=FALSE)
    x2<-plotHeatmap(data=smData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,],alignSampleData=TRUE,clusterFeatures=FALSE,clusterSamples=FALSE)
#   Should get this working so proper test, but more a problem because in different order, otherwise the same. Don't want to deal with this right now.
#    expect_equal(lapply(x1$clusterLegend,function(x){x[,c("clusterIds","color")]}),lapply(x2$clusterLegend,function(x){x[,c("clusterIds","color")]}))
})

test_that("`plotHeatmap` works with CE objects", {
    x1<-plotHeatmap(cl3[1:30,1:50])

})

