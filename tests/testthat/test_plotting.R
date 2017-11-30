context("Non-heatmap related plot functions")

source("create_objects.R")

test_that("`plotClusters` works with matrix, ClusterExperiment objects", {

    #test matrix version
    x<-plotClusters(object=clusterMatrix(ceSim))
    expect_equal(dim(clusterMatrix(ceSim)),dim(x$colors))
    expect_equal(dim(clusterMatrix(ceSim)),dim(x$aligned))
    expect_equal(length(x$clusterLegend),ncol(clusterMatrix(ceSim)))
    expect_error(plotClusters(object=clusterMatrix(ceSim),whichClusters="garbage"),"unable to find an inherited method")
    expect_error(plotClusters(object=clusterMatrix(ceSim),whichClusters=c(1,3,4)),"unable to find an inherited method")

    #test CE version
    x<-plotClusters(ceSim)
    expect_is(x,"ClusterExperiment")
    expect_equal( x,ceSim)

    xx<-plotClusters(ceSim,whichClusters="clusterMany")
    xx2<-plotClusters(object=ceSim,whichClusters="workflow") #only clusterMany values so should be the same
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
    x1<-plotClusters(object=ceSim,whichClusters="workflow",resetColors=TRUE)
    x2<-plotClusters(object=removeClusters(ceSim,"User"),resetColors=TRUE)
    whP<-clusterExperiment:::.TypeIntoIndices(ceSim,"workflow")
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
  x<-plotClusters(object=clusterMatrix(ceSim),sampleData=as.data.frame(colData(ceSim)))
  expect_equal(ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)),ncol(x$colors))
  expect_equal(ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)),ncol(x$aligned))
  expect_equal(length(x$clusterLegend),ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)))

  #test CE version
  test<- clusterMany(simCount,dimReduce="PCA",nDimReduce=c(5,10,50), isCount=TRUE,
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


test_that("plotClustersWorkflow", {
	cc<-clusterMany(mat, ks=c(3,4),nFilter=c(10,15),nPCADim=c(3,4),dimReduce=c("none","PCA","var"),clusterFunction="pam",
	                       subsample=FALSE, sequential=FALSE,run=TRUE,verbose=FALSE,
	                       isCount=FALSE)
	cc<-combineMany(cc,proportion=.7,whichClusters = "clusterMany")
	plotClustersWorkflow(cc)
	plotClustersWorkflow(cc,clusterManyLabels=FALSE)
	plotClustersWorkflow(cc,sortBy="clusterMany")
	plotClustersWorkflow(cc,sortBy="clusterMany",highlightOnTop=FALSE)
	plotClustersWorkflow(cc,highlightOnTop=FALSE)
	plotClustersWorkflow(cc,clusterManyLabels=FALSE,clusterLabels="test")
	expect_error(plotClustersWorkflow(cc,clusterManyLabels=c("1","2"),clusterLabels="test"),"number of cluster labels given in clusterManyLabels")
	expect_error(plotClustersWorkflow(cc,clusterManyLabels=TRUE,clusterLabels=c("A","test")),"number of cluster labels given in clusterLabels")

})


test_that("plotting helpers", {
  convertClusterLegend(smSimCE,output="aheatmap")
  convertClusterLegend(smSimCE,output="plotAndLegend")
  convertClusterLegend(smSimCE,output="matrixNames")
  convertClusterLegend(smSimCE,output="matrixColors")
  convertClusterLegend(smSimCE,output="matrixNames",whichClusters=c("cluster1"))
  convertClusterLegend(smSimCE,output="matrixNames",whichClusters=1:3)
  convertClusterLegend(smSimCE,output="plotAndLegend",whichClusters=c("cluster1"))
  expect_error(convertClusterLegend(smSimCE,output="plotAndLegend",whichClusters=1:3),"given whichClusters indicates more than 1 clustering which is not allowed for option")
  
  plotClusterLegend(smSimCE)
  plotClusterLegend(smSimCE,whichCluster="cluster1")
    showPalette()
  showPalette(massivePalette)
})


test_that("`plotBarplot` works with matrix, ClusterExperiment objects", {

    #test numeric matrix version
    plotBarplot(object=clusterMatrix(ceSim)[,1:2])
    #test vector version
    plotBarplot(object=clusterMatrix(ceSim)[,1])
    #check error
    expect_error(plotBarplot(object=clusterMatrix(ceSim)),"if 'object' a matrix, must contain at most 2 clusters")
    
    #test CE version with no defaults
    plotBarplot(ceSim)
    #test CE version whichClusters arguments
    plotBarplot(ceSim,whichClusters="workflow")
    plotBarplot(ceSim,whichClusters="primaryCluster")
    plotBarplot(ceSim)

    
    test<-ceSim
    clusterLegend(test)[[1]][,"name"]<-LETTERS[1:nrow(clusterLegend(ceSim)[[1]])]
    #test character matrix version
    plotBarplot(object=convertClusterLegend(test,output="matrixNames")[,1:2])
    #test character vector version
    plotBarplot(object=convertClusterLegend(test,output="matrixNames")[,1])
    #test labels argument
    plotBarplot(test,whichClusters=1:2,labels="id")
    plotBarplot(test,whichClusters=1:2,labels="name")
    #plotBarplot(ceSim,whichClusters="primaryCluster")
    
})

test_that("plotDimReduce works",{
	expect_silent(cl <- clusterMany(simData, nDimReduce=c(5, 10, 50), dimReduce="PCA",
	clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
	removeSil=c(TRUE,FALSE)))
	expect_silent(plotDimReduce(cl,legend="bottomright"))
	expect_silent(plotDimReduce(cl,legend=TRUE))
	expect_silent(clusterLegend(cl)[["nDimReduce=10,k=4,findBestK=FALSE,removeSil=TRUE"]][,"name"]<-LETTERS[1:5])
	expect_silent(plotDimReduce(cl,whichCluster="nDimReduce=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE))
	
	#test on object that doesn't have saved:
	expect_silent(clD<-plotDimReduce(ceSimData,dimReduce="PCA"))
	expect_equal(NCOL(reducedDim(clD,type="PCA")),2) #default.
	
	#higher dims.
	expect_silent(plotDimReduce(cl,whichCluster="nDimReduce=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE,whichDims=1:4))
	expect_error(plotDimReduce(cl,whichCluster="nDimReduce=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE,whichDims=158:200),"Invalid value for whichDims: larger than row or column")
	#force it to recalculate:
	expect_silent(plotDimReduce(cl,whichCluster="nDimReduce=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE,whichDims=51:58))
	
	
	
})