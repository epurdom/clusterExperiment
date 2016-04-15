library(clusterExperiment)
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version. Can delete this once package is stabilized.

test_that("`mergeClusters` works with matrix and ClusterExperiment objects", {
  clustNothing <- clusterMany(simData, ks=c(9,11),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE)
  clustCombined <- combineMany(clustNothing, whichClusters = "clusterMany")

  clustWithDendro <- makeDendrogram(clustCombined)
  mergedList <- mergeClusters(x=transform(clustCombined), countData=FALSE,
                              cl=primaryCluster(clustCombined),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", plotType="mergeMethod")

  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP")
  clustMerged
  })
