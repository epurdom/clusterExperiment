context("mergeCLusters")
source("create_objects.R")

test_that("`mergeClusters` works with matrix and ClusterExperiment objects", {
  cl1 <- clusterSingle(smSimData, clusterFunction="pam",
                       subsample=FALSE, sequential=FALSE,
                       clusterDArgs=list(k=6),isCount=FALSE)
  clustWithDendro <- makeDendrogram(cl1)
  #matrix version
  mergedList <- mergeClusters(x=transform(cl1), countData=FALSE,
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", plotType="mergeMethod")

  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none",plotType="all")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="locfdr", plotType="mergeMethod")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="MB", plotType="mergeMethod")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="JC", plotType="mergeMethod")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotType="mergeMethod")

  expect_true("mergeClusters" %in% clusterType(clustMerged))
  expect_true("mergeClusters" %in% colnames(clusterMatrix(clustMerged)))

  expect_warning(mergeClusters(x=transform(clustWithDendro), countData=FALSE,
                               cl=primaryCluster(clustWithDendro),plot="none",
                               mergeMethod="adjP",
                               dendro=clustWithDendro@dendro_samples),
                 "not equal to the number")

  #test if already exists
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP")
  primaryClusterIndex(clustMerged)<-2
  clustMerged<- makeDendrogram(clustMerged)
  clustMerged2<-mergeClusters(clustMerged,mergeMethod="adjP")
  expect_true("mergeClusters_1" %in% clusterType(clustMerged2))
  expect_true(!"combineMany_1" %in% clusterType(clustMerged2))
  expect_true(!"clusterMany_1" %in% clusterType(clustMerged2))
})
