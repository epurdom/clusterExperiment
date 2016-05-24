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
    
    bigCE<-ceSim
    bigCE<-makeDendrogram(bigCE,whichCluster="cluster1")
    x1<-bigCE
    #--- check clusterMany updates dendrogram correctly
    bigCE<-clusterMany(bigCE,k=2:8,clusterFunction="hierarchicalK")
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    #takes a long time!
    #expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    expect_error(makeDendrogram(bigCE,whichCluster="workflow"),"'whichCluster' must identify only a single clustering")
 
    #--- check combineMany updates dendrogram correctly
    bigCE<-combineMany(bigCE,proportion=0.3)
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    #expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    makeDendrogram(bigCE,whichCluster="combineMany") 
    
    
    #--- check mergeClusters updates dendrogram correctly
    bigCE<-mergeClusters(bigCE,mergeMethod="adjP",cutoff=0.2)
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    #expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    
    expect_error(getBestFeatures(bigCE,contrastType="Dendro"),"Primary cluster does not match the cluster on which the dendrogram was made")
})
test_that("`makeDendrogram` with dimReduce options", {
    x<-makeDendrogram(ccSE,dimReduce="PCA",ndims=3)
    expect_error(makeDendrogram(ccSE,dimReduce=c("PCA","var"),ndims=3))
    x2<-makeDendrogram(ccSE,dimReduce=c("PCA"),ndims=3,ignoreUnassigned=TRUE)
    expect_equal(x,x2)
    makeDendrogram(ccSE,dimReduce=c("var"),ndims=3,ignoreUnassigned=FALSE)
    makeDendrogram(ccSE,dimReduce=c("var"),ndims=3,ignoreUnassigned=TRUE)
    makeDendrogram(ccSE,dimReduce=c("cv"),ndims=3,ignoreUnassigned=FALSE)
    makeDendrogram(ccSE,dimReduce=c("cv"),ndims=3,ignoreUnassigned=TRUE)
    makeDendrogram(ccSE,dimReduce=c("mad"),ndims=3,ignoreUnassigned=FALSE)
    makeDendrogram(ccSE,dimReduce=c("mad"),ndims=3,ignoreUnassigned=TRUE)
    
})
