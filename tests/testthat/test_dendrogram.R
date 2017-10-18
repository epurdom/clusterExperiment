context("Dendrogram")
# library(devtools)
# load_all()
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
    
    expect_error(getBestFeatures(bigCE,contrastType="Dendro"),"only single cluster in clustering -- cannot run getBestFeatures")
	primaryClusterIndex(bigCE)<-3
	expect_error( getBestFeatures(bigCE,contrastType="Dendro"),"Primary cluster does not match the cluster on which the dendrogram was made")
})

test_that("plotDendrogram works with outgroup", {
    leg<-clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]
    leg[,"name"]<-letters[1:nrow(leg)]
    clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]<-leg
  dend <- makeDendrogram(ccSE)
  plotDendrogram(dend)
  plotDendrogram(dend,show.node.label=TRUE)
  plotDendrogram(dend,leafType="samples",labelType="name")
  plotDendrogram(dend,leafType="samples",labelType="name",removeOutbranch=FALSE)
  plotDendrogram(dend,leafType="samples",labelType="colorblock")
  plotDendrogram(dend,leafType="clusters",labelType="colorblock")
  plotDendrogram(dend,leafType="clusters",labelType="name")
  
  ## make all -2
  cl<-clusterMatrix(ccSE)[,1]
  cl[1]<- -2
  dend2<-addClusters(ccSE,cl,clusterLabel="newCluster")
  primaryClusterIndex(dend2)<-3
  dend2 <- makeDendrogram(dend2)
  plotDendrogram(dend2,leafType="clusters",labelType="colorblock")
  plotDendrogram(dend2,leafType="samples",labelType="colorblock")
  plotDendrogram(dend2,leafType="samples",labelType="colorblock",removeOutbranch=FALSE)

  ## make only single sample -2
  cl<-clusterMatrix(ccSE)[,1]
  cl[1]<-1
  dend3<-addClusters(ccSE,cl,clusterLabel="newCluster")
  primaryClusterIndex(dend3)<-3
  dend3 <- makeDendrogram(dend3)
  plotDendrogram(dend3,leafType="clusters",labelType="colorblock")
  plotDendrogram(dend3,leafType="samples",labelType="colorblock")
  plotDendrogram(dend3,leafType="samples",labelType="colorblock",removeOutbranch=FALSE)

  # This test breaks something. Needs to be figured out. 
  # ## make all -1 but two samples
  # ## can't be only 1 sample because then only 1 cluster so can't make a dendrogram...
  # cl<-rep(-1,length=nSamples(ccSE))
  # cl[1]<-3
  # cl[2]<-1
  # dend4<-addClusters(ccSE,cl,clusterLabel="missingCluster")
  # primaryClusterIndex(dend4)<-3
  # dend4 <- makeDendrogram(dend4)
  # plotDendrogram(dend4,leafType="clusters",labelType="colorblock")
  # plotDendrogram(dend4,leafType="samples",labelType="colorblock")
  # plotDendrogram(dend4,leafType="samples",labelType="colorblock",removeOutbranch=FALSE)

  ## make all -1 but one sample -- should get error bc only 1 cluster, can't make dendrogram; 
  ## in case this changes, this test will catch that need to fix plotDendrogram, which makes assumption that not possible.
  cl<-rep(-1,length=nSamples(ccSE))
  cl[1]<-3
  dend5<-addClusters(ccSE,cl,clusterLabel="missingCluster")
  primaryClusterIndex(dend5)<-3
  expect_error(makeDendrogram(dend5,dimReduce="none"),"Only 1 cluster given. Can not make a dendrogram.")
  expect_error(plotDendrogram(dend5,leafType="clusters",labelType="colorblock"),"No dendrogram is found for this ClusterExperiment Object. Run makeDendrogram first.")


    
})

test_that("plotDendrogram works with whichClusters", {
    leg<-clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]
    leg[,"name"]<-letters[1:nrow(leg)]
    clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]<-leg
  dend <- makeDendrogram(ccSE)
  dend<-mergeClusters(dend)
  plotDendrogram(dend,whichClusters="all",leafType="samples",labelType="colorblock")
  
  
})


test_that("plotDendrogram works with cluster missing", {
    leg<-clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]
    leg[,"name"]<-letters[1:nrow(leg)]
    clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]<-leg
  dend <- makeDendrogram(ccSE,unassignedSamples = c("cluster"))
  plotDendrogram(dend)
  plotDendrogram(dend,show.node.label=TRUE)
  plotDendrogram(dend,leafType="samples",labelType="name")
  plotDendrogram(dend,leafType="samples",labelType="colorblock")
  plotDendrogram(dend,leafType="clusters",labelType="colorblock")
  plotDendrogram(dend,leafType="clusters",labelType="name")
  
  ## make all -2
  dend2<-dend
  mat<-clusterMatrix(dend2)
  mat[1,1]<- -2
  dend2@clusterMatrix<-mat
  leg<-dend2@clusterLegend[[1]]
  leg<-leg[-which(leg[,"clusterIds"]== -1),]
  dend2@clusterLegend[[1]]<-leg
  dend2 <- makeDendrogram(dend2,unassignedSamples = c("cluster"))
  plotDendrogram(dend2,leafType="clusters",labelType="colorblock")
  plotDendrogram(dend2,leafType="samples",labelType="colorblock")
  
})