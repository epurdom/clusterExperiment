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

test_that("`makeDendrogram` works with whichCluster", {
    x1<-makeDendrogram(ccSE,whichCluster="Cluster2")
    x2<-makeDendrogram(ccSE,whichCluster=2)
    expect_equal(x1,x2)
    
    bigCE<-ccSE
    bigCE<-makeDendrogram(bigCE,whichCluster="Cluster2")
    
    #check clusterMany updates dendrogram correctly
    bigCE<-clusterMany(bigCE,k=2:8,clusterFunction="hierarchicalK")
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    expect_error(makeDendrogram(bigCE,whichCluster="workflow"),"'whichCluster' must identify only a single clustering")
 
    #check combineMany updates dendrogram correctly
    
    
    #check mergeClusters updates dendrogram correctly
    
})
