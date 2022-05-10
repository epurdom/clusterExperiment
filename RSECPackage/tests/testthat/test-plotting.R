context("Test Plotting from ClusterExperiment")

test_that("`plotClusters` works with RSECClass objects", {
  ## Check workflow
  expect_silent(xx<-plotClusters(ceSimCount,
    whichClusters="clusterMany"))
  expect_silent(xx2<-plotClusters(object=ceSimCount,
    whichClusters="workflow")) #only clusterMany values so should be the same
  expect_equal(xx2,xx)

  par(mfrow=c(1,2)) #so can visually check if desired.
  expect_silent(xx3<-plotClusters(ceSimCount,
    resetOrderSamples=TRUE,resetColors=TRUE,resetNames=TRUE))
  plotClusters(xx3,existingColors="all")
  expect_false(isTRUE(all.equal(xx2,xx3))) #not a great test. Doesn't really say whether does it right, just whether it does something!
  
  
  #RSEC object with mixture of workflow and other types
  x1<-plotClusters(object=ceSimCount,
    whichClusters="workflow",resetColors=TRUE)
  x2<-plotClusters(object=removeClusterings(ceSimCount,"User"),
    resetColors=TRUE)
  whP<-getClusterIndex(ceSimCount,"workflow")
  expect_equal(clusterLegend(x2),clusterLegend(x1)[whP])
  
  #no colData in test
  expect_silent( test<- clusterMany(simCount,
      reduceMethod="PCA", verbose=FALSE,
      nReducedDims=c(5,10,50), isCount=TRUE, makeMissingDiss=TRUE,
      clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE)
  ) )
  expect_error(plotClusters(test,
      colData=as.data.frame(colData(ceSimCount))),
      "no colData for object data")
  
})


test_that("plotClustersWorkflow", {
	expect_silent(cc<-clusterMany(mat, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),reduceMethod=c("none","PCA","var"),clusterFunction="pam",
	                       subsample=FALSE, sequential=FALSE,run=TRUE,verbose=FALSE,
	                       isCount=FALSE))
	expect_silent(cc<-makeConsensus(cc,proportion=.7,whichClusters = "clusterMany"))
	expect_silent(plotClustersWorkflow(cc))
	expect_silent(plotClustersWorkflow(cc,clusterManyLabels=FALSE))
	expect_silent(plotClustersWorkflow(cc,sortBy="clusterMany"))
	expect_silent(plotClustersWorkflow(cc,sortBy="clusterMany",highlightOnTop=FALSE))
	expect_silent(plotClustersWorkflow(cc,highlightOnTop=FALSE))
	expect_silent(plotClustersWorkflow(cc,clusterManyLabels=FALSE,clusterLabels="test"))
	expect_error(plotClustersWorkflow(cc,clusterManyLabels=c("1","2"),clusterLabels="test"),"number of cluster labels given in clusterManyLabels")
	expect_error(plotClustersWorkflow(cc,clusterManyLabels=TRUE,clusterLabels=c("A","test")),"number of cluster labels given in clusterLabels")

})
