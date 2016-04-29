context("mergeCLusters")
source("create_objects.R")

test_that("`mergeClusters` works with matrix and ClusterExperiment objects", {
  cl1 <- clusterSingle(smSimData, clusterFunction="pam",
                       subsample=FALSE, sequential=FALSE,
                       clusterDArgs=list(k=6),isCount=FALSE)
  clustWithDendro <- makeDendrogram(cl1)
  #matrix version
  mergedList <- mergeClusters(x=transform(cl1), isCount=FALSE,
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", plotType="mergeMethod")

  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none",plotType="all")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none", plotType="adjP")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none", plotType="locfdr")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="locfdr", plotType="mergeMethod")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="MB", plotType="mergeMethod")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="JC", plotType="mergeMethod")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotType="mergeMethod")
  expect_error(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none", plotType="mergeMethod"),"can only plot merge method values if one method is selected")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotType="none")
  
  expect_true("mergeClusters" %in% clusterType(clustMerged))
  expect_true("mergeClusters" %in% colnames(clusterMatrix(clustMerged)))

  expect_warning(mergeClusters(x=transform(clustWithDendro), isCount=FALSE,
                               cl=primaryCluster(clustWithDendro),plot="none",
                               mergeMethod="adjP",
                               dendro=clustWithDendro@dendro_samples),
                 "not equal to the number")

  #test if already exists
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP")
  primaryClusterIndex(clustMerged)<-2
  clustMerged<- makeDendrogram(clustMerged)
  clustMerged2<-mergeClusters(clustMerged,mergeMethod="adjP")
  expect_true("mergeClusters.1" %in% clusterType(clustMerged2))
  expect_true(!"combineMany.1" %in% clusterType(clustMerged2))
  expect_true(!"clusterMany.1" %in% clusterType(clustMerged2))
})

test_that("`mergeClusters` preserves the colData and rowData of SE", {

  cl <- clusterSingle(smSimSE, clusterFunction="pam",
                       subsample=FALSE, sequential=FALSE,
                       clusterDArgs=list(k=6),isCount=FALSE)
  cl <- makeDendrogram(cl)
  cl <- mergeClusters(cl, mergeMethod = "adjP")
  expect_equal(colData(cl),colData(smSimSE))
  expect_equal(rownames(cl),rownames(smSimSE))
  expect_equal(colnames(cl),colnames(smSimSE))
  expect_equal(metadata(cl),metadata(smSimSE))
  expect_equal(rowData(cl),rowData(smSimSE))

})
