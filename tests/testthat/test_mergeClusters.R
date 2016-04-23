context("mergeCLusters")
source("create_objects.R")

test_that("`mergeClusters` works with matrix and ClusterExperiment objects", {
  clustNothing <- clusterMany(simData, ks=c(9,11),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE)
  clustCombined <- combineMany(clustNothing, whichClusters = "clusterMany")

  clustWithDendro <- makeDendrogram(clustCombined)
  #matrix version
  mergedList <- mergeClusters(x=transform(clustCombined), countData=FALSE,
                              cl=primaryCluster(clustCombined),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", plotType="mergeMethod")

  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP")
  expect_true("mergeClusters" %in% clusterType(clustMerged))
  expect_true("mergeClusters" %in% colnames(clusterMatrix(clustMerged)))

  expect_warning(mergeClusters(x=transform(clustWithDendro), countData=FALSE,
                               cl=primaryCluster(clustWithDendro),
                               dendro=clustWithDendro@dendro_samples),
                 "not equal to the number")

  mergedList <- mergeClusters(x=simCount, countData=TRUE,
                              cl=primaryCluster(clustWithDendro),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP")

  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="locfdr")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="MB")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="JC")

  #test if already exists
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP")
  primaryClusterIndex(clustMerged)<-2
  clustMerged<- makeDendrogram(clustMerged)
  clustMerged2<-mergeClusters(clustMerged)
  expect_true("mergeClusters_1" %in% clusterType(clustMerged2))
  expect_true(!"combineMany_1" %in% clusterType(clustMerged2))
  expect_true(!"clusterMany_1" %in% clusterType(clustMerged2))
})
