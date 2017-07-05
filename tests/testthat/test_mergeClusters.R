context("mergeCLusters")
source("create_objects.R")

test_that("`mergeClusters` works with matrix and ClusterExperiment objects", {
  cl1 <- clusterSingle(smSimData, 
                       subsample=FALSE, sequential=FALSE,
                       mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),isCount=FALSE)
  leg<-clusterLegend(cl1)[[primaryClusterIndex(cl1)]]
  leg[,"name"]<-letters[1:6]
  clusterLegend(cl1)[[primaryClusterIndex(cl1)]]<-leg
  clustWithDendro <- makeDendrogram(cl1)
  #matrix version
  mergedList <- mergeClusters(x=transform(cl1), isCount=FALSE,
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", plotInfo="mergeMethod")

	#check plotting types:
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none",plotInfo="all")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none", plotInfo="adjP")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none", plotInfo="locfdr")
  expect_error(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none", plotInfo="mergeMethod"),"can only plot 'mergeMethod' results if one method is selected")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="none")

  #check all methods run
  for(method in clusterExperiment:::.availMergeMethods){
	  clustMerged <- mergeClusters(clustWithDendro, mergeMethod=method, plotInfo="mergeMethod")
  }
  
  expect_true("mergeClusters" %in% clusterTypes(clustMerged))
  expect_true("mergeClusters" %in% colnames(clusterMatrix(clustMerged)))

	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",labelType="colorblock")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",labelType="name")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",labelType="colorblock")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",labelType="name")

  expect_error(mergeClusters(x=transform(clustWithDendro), isCount=FALSE,
                               cl=primaryCluster(clustWithDendro),plot="none",
                               mergeMethod="adjP",
                               dendro=clustWithDendro@dendro_samples),
                 "Not a valid input dendrogram")

  #test if already exists
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP")
  primaryClusterIndex(clustMerged)<-2
  clustMerged<- makeDendrogram(clustMerged)
  clustMerged2<-mergeClusters(clustMerged,mergeMethod="adjP")
  expect_true("mergeClusters.1" %in% clusterTypes(clustMerged2))
  expect_true(!"combineMany.1" %in% clusterTypes(clustMerged2))
  expect_true(!"clusterMany.1" %in% clusterTypes(clustMerged2))
  removeClusters(clustMerged, whichRemove = "mergeClusters")
})

test_that("`mergeClusters` preserves the colData and rowData of SE", {

  cl <- clusterSingle(smSimSE, 
                       subsample=FALSE, sequential=FALSE,
                       mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),isCount=FALSE)
  cl <- makeDendrogram(cl)
  cl <- mergeClusters(cl, mergeMethod = "adjP")
  expect_equal(colData(cl),colData(smSimSE))
  expect_equal(rownames(cl),rownames(smSimSE))
  expect_equal(colnames(cl),colnames(smSimSE))
  expect_equal(metadata(cl),metadata(smSimSE))
  expect_equal(rowData(cl),rowData(smSimSE))

})


test_that("`mergeClusters` works with unassignedSamples", {

  clustWithDendro <- makeDendrogram(ceSim,unassignedSamples = c("outgroup"))

	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",labelType="colorblock")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",labelType="name")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",labelType="colorblock")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",labelType="name")

  clustWithDendro <- makeDendrogram(ceSim,unassignedSamples = c("cluster"))

	expect_warning(mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",labelType="colorblock"))
	expect_warning(mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",labelType="name"))
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",labelType="colorblock")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",labelType="name")


})

