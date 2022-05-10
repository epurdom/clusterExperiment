context("Test Plotting from ClusterExperiment")

test_that("`plotClusters` works with RSECClass objects", {
  ## Check workflow
  xx<-plotClusters(ceSimCount,whichClusters="clusterMany")
  xx2<-plotClusters(object=ceSimCount,whichClusters="workflow") #only clusterMany values so should be the same
  expect_equal(xx2,xx)
  
  
  #RSEC object with mixture of workflow and other types
  x1<-plotClusters(object=ceSimCount,
    whichClusters="workflow",resetColors=TRUE)
  x2<-plotClusters(object=removeClusterings(ceSimCount,"User"),
    resetColors=TRUE)
  whP<-getClusterIndex(ceSimCount,"workflow")
  expect_equal(clusterLegend(x2),clusterLegend(x1)[whP])
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
