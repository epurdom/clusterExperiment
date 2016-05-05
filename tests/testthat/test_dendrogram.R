context("Dendrogram")
source("create_objects.R")

test_that("`makeDendrogram` works with matrix, ClusterExperiment objects", {
    #test matrix version
    makeDendrogram(mat, primaryCluster(cc))
    makeDendrogram(mat, primaryCluster(cc), unassigned="cluster")
    makeDendrogram(mat, primaryCluster(cc), unassigned="remove")

    #test CE version
    makeDendrogram(cc)
    makeDendrogram(cc, unassigned="cluster")
    expect_error(makeDendrogram(cc, unassigned="remove"))

    #test matrix version
    expect_equal(nobs(makeDendrogram(mat, primaryCluster(cc), unassigned="remove")$samples),
                 length(primaryCluster(cc))-2)

    #test proper error if only single cluster:
    fakeCluster<-rep(1,nSamples(cc))
    cc1<-addClusters(cc,fakeCluster)
    primaryClusterIndex(cc1)<-3
    expect_error(makeDendrogram(cc1),"Only 1 cluster given")
    fakeCluster[c(1:2)]<- -1
    cc1<-addClusters(cc,fakeCluster)
    primaryClusterIndex(cc1)<-3
    expect_error(makeDendrogram(cc1),"Only 1 cluster given")
})

test_that("`makeDendrogram` preserves the colData and rowData of SE", {

  dend <- makeDendrogram(ccSE)

  expect_equal(colData(dend),colData(se))
  expect_equal(rownames(dend),rownames(se))
  expect_equal(colnames(dend),colnames(se))
  expect_equal(metadata(dend),metadata(se))
  expect_equal(rowData(dend),rowData(se))

})
