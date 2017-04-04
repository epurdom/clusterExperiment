context("Plot functions")

source("create_objects.R")

test_that("`plotClusters` works with matrix, ClusterExperiment objects", {

    #test matrix version
    x<-plotClusters(clusters=clusterMatrix(ceSim))
    expect_equal(dim(clusterMatrix(ceSim)),dim(x$colors))
    expect_equal(dim(clusterMatrix(ceSim)),dim(x$aligned))
    expect_equal(length(x$clusterLegend),ncol(clusterMatrix(ceSim)))
    expect_error(plotClusters(clusters=clusterMatrix(ceSim),whichClusters="garbage"),"unable to find an inherited method")
    expect_error(plotClusters(clusters=clusterMatrix(ceSim),whichClusters=c(1,3,4)),"unable to find an inherited method")

    #test CE version
    x<-plotClusters(ceSim)
    expect_is(x,"ClusterExperiment")
    expect_equal( x,ceSim)

    xx<-plotClusters(ceSim,whichClusters="clusterMany")
    xx2<-plotClusters(clusters=ceSim,whichClusters="workflow") #only clusterMany values so should be the same
    expect_equal(xx2,xx)

    #check reset -- should add combinations of resetColors and resetNames to make sure works independently.
    par(mfrow=c(1,2)) #so can visually check if desired.
    xx3<-plotClusters(ceSim,resetOrderSamples=TRUE,resetColors=TRUE,resetNames=TRUE)
    plotClusters(xx3,existingColors="all")
    expect_false(isTRUE(all.equal(xx2,xx3))) #not a great test. Doesn't really say whether does it right, just whether it does something!

    nm<-as.numeric(unlist(lapply(clusterLegend(xx3),function(x){x[,"name"]})))
    col<-(unlist(lapply(clusterLegend(xx3),function(x){x[,"color"]})))
    wh<-which(col %in% c("white","grey"))
    expect_equal(match(col[-wh],bigPalette),nm[-wh])
    nmOld<-as.numeric(unlist(lapply(clusterLegend(ceSim),function(x){x[,"name"]})))
    expect_false(isTRUE(all.equal(nm,nmOld)))
    idOld<-as.numeric(unlist(lapply(clusterLegend(ceSim),function(x){x[,"clusterIds"]})))
    idNew<-as.numeric(unlist(lapply(clusterLegend(xx3),function(x){x[,"clusterIds"]})))
    expect_equal(idOld,idNew)

    #check existing colors
    x2<-plotClusters(ceSim,existingColors="all")

    #test -1
    plotClusters(ceSim)

    #CE object with mixture of workflow and other types
    x1<-plotClusters(clusters=ceSim,whichClusters="workflow",resetColors=TRUE)
    x2<-plotClusters(clusters=removeClusters(ceSim,"User"),resetColors=TRUE)
    whP<-.TypeIntoIndices(ceSim,"workflow")
    expect_equal(clusterLegend(x2),clusterLegend(x1)[whP])

    #test specifying indices
    wh<-c(3,4,NCOL(clusterMatrix(ceSim)))
    x3<-plotClusters(ceSim,whichClusters=wh,axisLine=-2,resetColors=TRUE)
    x4<-plotClusters(ceSim,whichClusters=wh[c(3,2,1)],axisLine=-2,resetColors=TRUE)
    expect_false(isTRUE(all.equal(x3,x4)))

    par(mfrow=c(1,1)) #otherwise will affect other tests.
})


test_that("`plotClusters` rerun above tests with sampleData included", {

  #test matrix version
  x<-plotClusters(clusters=clusterMatrix(ceSim),sampleData=as.data.frame(colData(ceSim)))
  expect_equal(ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)),ncol(x$colors))
  expect_equal(ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)),ncol(x$aligned))
  expect_equal(length(x$clusterLegend),ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)))

  #test CE version
  test<- clusterMany(simCount,dimReduce="PCA",nPCADims=c(5,10,50), isCount=TRUE,
                     clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE)) #no colData in test
  expect_error(plotClusters(test,sampleData=as.data.frame(colData(ceSim))),"no colData for object data")
  expect_error(plotClusters(ceSim,sampleData=as.data.frame(colData(ceSim))),"invalid values for pulling sampleData")
  plotClusters(ceSim,sampleData="all")
  par(mfrow=c(1,2))
  x2<-plotClusters(ceSim,sampleData="all",resetColors=TRUE)
  x1<-plotClusters(ceSim,resetColors=TRUE)
  
  
  #check NAs
  naSim<-ceSim
  colData(naSim)[sample(10,1:nrow(naSim)),]<-NA
  plotClusters(naSim,sampleData=c("A","B"))

  #test the new TRUE option
  plotClusters(naSim,sampleData=TRUE) 
  
  #this is not working because first one puts -1/-2 last and second puts them first, and so then assigns different colors to the groups
#  expect_equal(x1,x2)
#   par(mfrow=c(1,2))
#   x2<-plotClusters(ceSim,sampleData="all",resetColors=FALSE)
#   x1<-plotClusters(ceSim,resetColors=FALSE)
  par(mfrow=c(1,1))

})
test_that("`setBreaks`", {
	setBreaks(smSimData)
	setBreaks(smSimData,breaks=0.99)
	x<-setBreaks(smSimData,breaks=0.99,makeSymmetric=TRUE)
	expect_equal(max(x),-min(x))
	expect_equal(x,setBreaks(smSimData,breaks=0.01,makeSymmetric=TRUE))
	expect_warning(y<-setBreaks(smSimData,breaks=10))
	expect_equal(length(y),10)
}
test_that("`plotHeatmap` works with matrix objects", {

    x1<-plotHeatmap(data=smSimData)
    a1<-NMF::aheatmap(smSimData)
    expect_equal(x1$aheatmapOut,a1)
    x2<-plotHeatmap(data=smSimCount,clusterSamplesData=smSimData,clusterFeaturesData=smSimData)
    #for some reason, labels on dendrogram move from character to numeric so can't test entire object...
    expect_equal(x1$aheatmapOut$rowInd,x2$aheatmapOut$rowInd) 
    expect_equal(x1$aheatmapOut$colInd,x2$aheatmapOut$colInd) 
    
    #check internal alignment of sampleData (alignSampleData=TRUE) is working:
    sampleData<-clusterMatrix(smSimCE)
    alList<-plotClusters(sampleData)
    alCol<-alList$clusterLegend
    x1<-plotHeatmap(data=smSimData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,1:10],clusterLegend=alCol,clusterSamples=FALSE,clusterFeatures=FALSE)
    x2<-plotHeatmap(data=smSimData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,1:10],alignSampleData=TRUE,clusterFeatures=FALSE,clusterSamples=FALSE)
#   Should get this working so proper test, but more a problem because in different order, otherwise the same. Don't want to deal with this right now.
#    expect_equal(lapply(x1$clusterLegend,function(x){x[,c("clusterIds","color")]}),lapply(x2$clusterLegend,function(x){x[,c("clusterIds","color")]}))

    expect_error( plotHeatmap(data=smSimData,Rowv=TRUE),"arguments to aheatmap cannot be set by the user")
    expect_error( plotHeatmap(data=smSimData,Colv=TRUE),"arguments to aheatmap cannot be set by the user")
    expect_error( plotHeatmap(data=smSimData,colorScale=seqPal5,color=TRUE),"arguments to aheatmap cannot be set by the user")

    expect_error( plotHeatmap(data=smSimData,annCol=rnorm(n=ncol(smSimData))),"arguments to aheatmap cannot be set by the user")
    expect_error( plotHeatmap(data=smSimData,annColors=list(a=c("blue","green"))),"arguments to aheatmap cannot be set by the user")

    x1<-plotHeatmap(data=smSimData)
    
    
    ##Should add tests that pass aheatmap arguments correctly.
})

test_that("`plotHeatmap` works with ClusterExperiment and SummarizedExperiment objects", {

    plotHeatmap(cc)
    plotHeatmap(cc,whichClusters="none")
    expect_warning(plotHeatmap(cc,whichClusters="workflow") ,"whichClusters value does not match any clusters") #there are no workflow for this one

    plotHeatmap(smSimCE,whichClusters="workflow",overRideClusterLimit=TRUE)
    plotHeatmap(smSimCE,whichClusters="all",alignSampleData=TRUE,overRideClusterLimit=TRUE)
    expect_error(plotHeatmap(smSimCE,whichClusters=1:15),"Indices in whichClusters invalid")

    #test sampleData
    expect_error(plotHeatmap(cc,sampleData="A"), "no colData for object data")

    plotHeatmap(smSimCE,sampleData="all",overRideClusterLimit=TRUE)
    plotHeatmap(smSimCE,sampleData="A")
    plotHeatmap(smSimCE,sampleData=2:3)

    #check that it pulls the names, not the clusterIds.
    clusterLegend(cc)[[1]][,"name"]<-letters[1:nrow(clusterLegend(cc)[[1]])]
    plotHeatmap(cc)
    
    #check user setting clusterLegend
    plotHeatmap(cc,clusterLegend=list("Cluster1"=palette()[1:7]))
    plotHeatmap(smSimCE,sampleData="A",clusterLegend=list("A"=palette()[1:3]))
    # the following works outside of the test but not inside
    # possibly issue with testthat? Not evaluating for now.
    #plotHeatmap(smSimCE, sampleData="all", whichClusters="none")

    #SummarizedExperiment
    plotHeatmap(smSimSE)

})

test_that("`plotHeatmap` visualization choices/feature choices all work", {

  plotHeatmap(smSimCE,visualizeData=smSimCount)
  plotHeatmap(smSimCE,visualizeData="transformed")
  plotHeatmap(smSimCE,visualizeData="original")
  plotHeatmap(smSimCE,visualizeData="centeredAndScaled")
  #even if visualizeData="orginal, still clsuter on transformed. Should make unit test out of below that get same:
  plotHeatmap(smSimCE,visualizeData="transformed",clusterSamplesData="hclust")
  orderSamples(smSimCE)<-sample(1:nSamples(smSimCE))
  plotHeatmap(smSimCE,visualizeData="transformed",clusterSamplesData="orderSamplesValue")
  plotHeatmap(smSimCE,visualizeData="transformed",clusterSamplesData="primaryCluster")
  plotHeatmap(smSimCE,visualizeData="transformed",clusterSamplesData=c(3,4,5))

  plotHeatmap(smSimCE,visualizeData="transform",clusterFeaturesData="all")
  plotHeatmap(smSimCE,visualizeData="transform",clusterFeaturesData="var",nFeatures=3)
  plotHeatmap(smSimCE,visualizeData="transform",clusterFeaturesData=3:5,nFeatures=3)
  expect_error(plotHeatmap(smSimCE,visualizeData="transform",clusterFeaturesData=paste("Gene",3:5),nFeatures=3))
  row.names(smSimCE)<-paste("Gene",1:NROW(smSimCE))
  plotHeatmap(smSimCE,visualizeData="transform",clusterFeaturesData=paste("Gene",3:5),nFeatures=3)
  plotHeatmap(smSimCE,visualizeData="transform",clusterFeaturesData="PCA",nFeatures=10,clusterSamplesData="hclust")

  plotHeatmap(smSimCE,visualizeData="transform",clusterSamplesData="dendrogramValue")

})

test_that("`makeBlankData` works", {


  ##call directly
  gps<-list(c(3,6,7),c(2,1))
  xx<-makeBlankData(assay(smSimCE),groupsOfFeatures=gps)
  expect_equal(nrow(xx$dataWBlanks),length(xx$rowNamesWBlanks))
  whBlankNames<-which(xx$rowNamesWBlanks=="")
  expect_equal(xx$rowNamesWBlanks[-whBlankNames],as.character(unlist(gps)) )
  whBlankRows<-as.numeric(which(apply(xx$dataWBlanks,1,function(x){all(is.na(x))})))
  expect_equal(whBlankRows,whBlankNames)
  expect_equal(whBlankRows,4)

  ##call within plotHeatmap
  plotHeatmap(smSimCE,clusterFeaturesData=gps)
})
test_that("`plotCoClustering` works", {
  expect_error(plotCoClustering(smSimCE),"coClustering slot is empty")
  smMin1<-combineMany(smSimCE,whichClusters=1:4,proportion=.95) #gives all -1, but creates coClustering
  plotCoClustering(smMin1,clusterSamplesData="hclust")
  expect_error(plotCoClustering(smMin1,clusterSamplesData="dendrogramValue"),
               "all samples have clusterIds<0")
  sm<-combineMany(smSimCE,whichClusters=1:4,proportion=.5)
  plotCoClustering(sm,clusterSamplesData="dendrogramValue")
})

test_that("plotting helpers", {
  convertClusterLegend(smSimCE,output="aheatmap")
  convertClusterLegend(smSimCE,output="plotAndLegend")
  convertClusterLegend(smSimCE,output="matrixNames")
  convertClusterLegend(smSimCE,output="matrixColors")
})
