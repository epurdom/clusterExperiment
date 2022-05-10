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