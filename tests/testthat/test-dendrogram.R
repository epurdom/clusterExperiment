context("Dendrogram")
# library(devtools)
# load_all()


test_that("`makeDendrogram` works with matrix",{
    #test matrix version
    expect_silent(makeDendrogram(mat, primaryCluster(cc)))
    expect_silent(makeDendrogram(mat, primaryCluster(cc), unassigned="cluster"))
    expect_silent(makeDendrogram(mat, primaryCluster(cc), unassigned="remove"))

 
    #test matrix version
    expect_equal(phylobase::nTips(makeDendrogram(mat, primaryCluster(cc), unassigned="remove")$samples),
                 length(primaryCluster(cc))-2)
})
			 
test_that("`makeDendrogram` works with ClusterExperiment objects", {
    #test CE version
    expect_silent(makeDendrogram(cc))
    expect_silent(makeDendrogram(cc, unassigned="cluster"))
    expect_error(makeDendrogram(cc, unassigned="remove"),"should be one of") #not valid option

    #test proper error if only single cluster:
    fakeCluster<-rep(1,nSamples(cc))
    expect_silent(cc1<-addClusterings(cc,fakeCluster))
    expect_silent(primaryClusterIndex(cc1)<-3)
    expect_error(makeDendrogram(cc1),"Only 1 cluster given")
    fakeCluster[c(1:2)]<- -1
    expect_silent(cc1<-addClusterings(cc,fakeCluster))
    expect_silent(primaryClusterIndex(cc1)<-3)
    expect_error(makeDendrogram(cc1),"Only 1 cluster given")
	
	#check whichClusters for clusterMatrix
 	expect_silent(dend<-makeDendrogram(ceSim))
	expect_equal(ncol(clusterMatrix(dend,whichClusters="dendro")),1)
	
	#--------
    #Check filters created when ignoriing cluster for variance
	#--------
    expect_silent(cl1 <- clusterSingle(smSimData, subsample=FALSE, sequential=FALSE, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),isCount=FALSE))
    leg<-clusterLegend(cl1)[[primaryClusterIndex(cl1)]]
    leg[,"name"]<-letters[1:6]
    expect_silent(clusterLegend(cl1)[[primaryClusterIndex(cl1)]]<-leg)
    expect_silent(clustWithDendro <- makeDendrogram(cl1))
	expect_equal(filterNames(clustWithDendro),"mad_clusterSingle")
    #redo
	expect_silent(clustWithDendro <- makeDendrogram(clustWithDendro))
	expect_equal(filterNames(clustWithDendro),"mad_clusterSingle")
    
	
	#--------
    #Check remove clusters when have dendrogram
	#--------
	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adj", plotInfo="none",DEMethod="limma"),"Merging will be done on ' clusterSingle ', with clustering index 1")
    expect_silent(removeClusterings(clustMerged,whichClusters="mergeClusters")) #remove merged, keep one with dendrogram
    expect_silent(removeClusterings(clustMerged,whichClusters=2)) #remove one with dendrogram
  
    #--------
	#test subsetting works if have dendrogram attached:
    #--------
	expect_silent(x<-makeDendrogram(ccSE,reduceMethod="PCA",nDims=3))
    expect_silent(x[1:3,1:2])
  
})

test_that("`makeDendrogram` works with hdf5",{
    expect_silent(clustNothing <- clusterMany(hdfSCE,
		 ks=c(3,4),clusterFunction="pam",
		 subsample=FALSE, sequential=FALSE,
		 isCount=FALSE,verbose=FALSE))
    expect_silent(clustNothing<-makeConsensus(clustNothing, proportion=1,whichClusters = "clusterMany"))
	expect_silent(makeDendrogram(clustNothing))
	
})

test_that("`makeDendrogram` preserves the colData and rowData of SE", {

  expect_silent(whCl<-primaryClusterIndex(ccSE))
  expect_silent(dend <- makeDendrogram(ccSE,reduceMethod="mad", whichCluster=whCl))

  expect_equal(colData(dend),colData(se))
  expect_equal(rownames(dend),rownames(se))
  expect_equal(colnames(dend),colnames(se))
  expect_equal(metadata(dend),metadata(se))
  x<-rowData(dend)
  expect_silent(newFilter<-clusterExperiment:::.makeClusterFilterStats("mad",clusterLabels(ccSE)[whCl]))
  expect_equal(x[,-grep(newFilter,colnames(x))],rowData(se))

})

test_that("`makeDendrogram` with reduceMethod options", {
    expect_silent(x<-makeDendrogram(ccSE,reduceMethod="PCA",nDims=3))
	expect_error(makeDendrogram(ccSE,reduceMethod=c("PCA","var"),nDims=3))
    expect_silent(x2<-makeDendrogram(ccSE,reduceMethod=c("PCA"),nDims=3,filterIgnoresUnassigned=TRUE))
    expect_equal(x,x2)
    expect_silent(makeDendrogram(ccSE,reduceMethod=c("var"),nDims=3,filterIgnoresUnassigned=FALSE))
    expect_silent(makeDendrogram(ccSE,reduceMethod=c("var"),nDims=3,filterIgnoresUnassigned=TRUE))
    expect_silent(makeDendrogram(ccSE,reduceMethod=c("abscv"),nDims=3,filterIgnoresUnassigned=FALSE))
    expect_silent(makeDendrogram(ccSE,reduceMethod=c("abscv"),nDims=3,filterIgnoresUnassigned=TRUE))
    expect_silent(makeDendrogram(ccSE,reduceMethod=c("mad"),nDims=3,filterIgnoresUnassigned=FALSE))
    expect_silent(makeDendrogram(ccSE,reduceMethod=c("mad"),nDims=3,filterIgnoresUnassigned=TRUE))
    
})
test_that("`makeDendrogram` works with whichCluster", {
    expect_silent(x1<-makeDendrogram(ccSE,whichCluster="Cluster2"))
    expect_silent(x2<-makeDendrogram(ccSE,whichCluster=2))
    expect_equal(x1,x2)
    
    bigCE<-ceSim
    expect_silent(bigCE<-makeDendrogram(bigCE,whichCluster="cluster1"))
    x1<-bigCE
    #--- check clusterMany updates dendrogram correctly
    expect_silent(bigCE<-clusterMany(bigCE,k=2:8,clusterFunction="hierarchicalK"))
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    #takes a long time!
    #expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    expect_error(makeDendrogram(bigCE,whichCluster="workflow"),"'whichCluster' must identify only a single clustering")
 
    #--- check makeConsensus updates dendrogram correctly
    expect_message(bigCE<-makeConsensus(bigCE,proportion=0.3),"no clusters specified to combine, using results from clusterMany")
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    expect_silent(makeDendrogram(bigCE,whichCluster="makeConsensus") )
    
    
    #--- check mergeClusters updates dendrogram correctly
    expect_message(bigCE<-mergeClusters(bigCE,mergeMethod="adjP",cutoff=0.2, DEMethod="limma"),"Merging will be done on ")
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],clusterLabels(x1)[x1@dendro_index])
    #Note, x1 doesn't give any merged clusters in dendrogram because before mergeClusters step...  so going to remove that element from bigCE
	tdf<-phylobase::tdata(bigCE@dendro_clusters)
	tdf$ClusterIdMerge<-NA
	phylobase::tdata(bigCE@dendro_clusters)<-tdf
	expect_equal(bigCE@dendro_clusters,x1@dendro_clusters)  
    expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    
    expect_error(getBestFeatures(bigCE,contrastType="Dendro"),"only single cluster in clustering -- cannot run getBestFeatures")
	expect_silent(primaryClusterIndex(bigCE)<-3)
	expect_error( getBestFeatures(bigCE,contrastType="Dendro"),"does not match either the cluster on which the dendrogram was made or the merge cluster from this dendrogram")
})

test_that("plotDendrogram works with colData", {
  leg<-clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]
  leg[,"name"]<-letters[1:nrow(leg)]
  clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]<-leg
	expect_silent(dend <- makeDendrogram(ccSE))
	expect_silent(plotDendrogram(dend,colData="A"))
	expect_warning(plotDendrogram(dend,colData=c("A","B","C")),"implies using columns of colData that are continuous")
	
	#note that legA doesn't give colors for everything -- only some. 
	legA<-leg[4:7,]
	legA[,"color"]<-tail(massivePalette,4)
	expect_silent(plotDendrogram(dend,colData="A",clusterLegend=list("A"=legA)))

	expect_silent(plotDendrogram(dend,colData=c("A","C"),clusterLegend=list("A"=legA)))
	
})

test_that("plotDendrogram works with outgroup", {
    leg<-clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]
    leg[,"name"]<-letters[1:nrow(leg)]
    clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]<-leg
  dend <- makeDendrogram(ccSE)
  expect_silent(plotDendrogram(dend))
  expect_silent(plotDendrogram(dend,show.node.label=TRUE))
  expect_silent(plotDendrogram(dend,leafType="samples",plotType="name"))
  expect_silent(plotDendrogram(dend, leafType="samples",plotType="name",removeOutbranch=FALSE))
  expect_silent(plotDendrogram(dend,leafType="samples",plotType="colorblock"))
  expect_silent(plotDendrogram(dend,leafType="clusters",plotType="colorblock"))
  expect_silent(plotDendrogram(dend,leafType="clusters",plotType="name"))
  
  ## make all -2
  cl<-clusterMatrix(ccSE)[,1]
  cl[1]<- -2
  dend2<-addClusterings(ccSE,cl,clusterLabel="newCluster")
  primaryClusterIndex(dend2)<-3
  dend2 <- makeDendrogram(dend2)
  expect_silent(plotDendrogram(dend2,leafType="clusters",plotType="colorblock"))
  expect_silent(plotDendrogram(dend2,leafType="samples",plotType="colorblock"))
  expect_silent(plotDendrogram(dend2,leafType="samples",plotType="colorblock",removeOutbranch=FALSE))

  ## make only single sample -2
  cl<-clusterMatrix(ccSE)[,1]
  cl[1]<-1
  expect_silent(dend3<-addClusterings(ccSE,cl,clusterLabel="newCluster"))
  expect_silent(primaryClusterIndex(dend3)<-3)
  expect_silent(dend3 <- makeDendrogram(dend3))
  expect_silent(plotDendrogram(dend3,leafType="clusters",plotType="colorblock"))
  expect_silent(plotDendrogram(dend3,leafType="samples",plotType="colorblock"))
  expect_silent(plotDendrogram(dend3,leafType="samples",plotType="colorblock",removeOutbranch=FALSE))

  # This test breaks something. Needs to be figured out. 
  # ## make all -1 but two samples
  # ## can't be only 1 sample because then only 1 cluster so can't make a dendrogram...
  # cl<-rep(-1,length=nSamples(ccSE))
  # cl[1]<-3
  # cl[2]<-1
  # dend4<-addClusterings(ccSE,cl,clusterLabel="missingCluster")
  # primaryClusterIndex(dend4)<-3
  # dend4 <- makeDendrogram(dend4)
  # plotDendrogram(dend4,leafType="clusters",plotType="colorblock")
  # plotDendrogram(dend4,leafType="samples",plotType="colorblock")
  # plotDendrogram(dend4,leafType="samples",plotType="colorblock",removeOutbranch=FALSE)

  ## make all -1 but one sample -- should get error bc only 1 cluster, can't make dendrogram; 
  ## in case this changes, this test will catch that need to fix plotDendrogram, which makes assumption that not possible.
  cl<-rep(-1,length=nSamples(ccSE))
  cl[1]<-3
  expect_silent(dend5<-addClusterings(ccSE,cl,clusterLabel="missingCluster"))
  expect_silent(primaryClusterIndex(dend5)<-3)
  expect_error(makeDendrogram(dend5,reduceMethod="none"),"Only 1 cluster given. Can not make a dendrogram.")
  expect_error(plotDendrogram(dend5,leafType="clusters",plotType="colorblock"),"No dendrogram is found for this ClusterExperiment Object. Run makeDendrogram first.")


    
})

test_that("plotDendrogram works with whichClusters", {
    leg<-clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]
    leg[,"name"]<-letters[1:nrow(leg)]
    clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]<-leg
  expect_silent(dend <- makeDendrogram(ccSE))
  expect_message(dend<-mergeClusters(dend,DEMethod="limma"))
  expect_silent(plotDendrogram(dend,whichClusters="all",leafType="samples",plotType="colorblock"))
  
  
})


test_that("plotDendrogram works with cluster missing", {
    expect_silent(leg<-clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]])
    leg[,"name"]<-letters[1:nrow(leg)]
    expect_silent(clusterLegend(ccSE)[[primaryClusterIndex(ccSE)]]<-leg)
  expect_silent(dend <- makeDendrogram(ccSE,unassignedSamples = c("cluster")))
  expect_silent(plotDendrogram(dend))
  expect_silent(plotDendrogram(dend,show.node.label=TRUE))
  expect_silent(plotDendrogram(dend,leafType="samples",plotType="name"))
  expect_silent(plotDendrogram(dend,leafType="samples",plotType="colorblock"))
  expect_silent(plotDendrogram(dend,leafType="clusters",plotType="colorblock"))
  expect_silent(plotDendrogram(dend,leafType="clusters",plotType="name"))
  
  ## make all the <0 ones -2 value
  dend2<-dend
  dmat<-clusterMatrix(dend2)
  dmat[1,1]<- -2
  dend2@clusterMatrix<-dmat
  leg<-dend2@clusterLegend[[1]]
  leg<-leg[-which(leg[,"clusterIds"]== -1),]
  dend2@clusterLegend[[1]]<-leg
  expect_silent(dend2 <- makeDendrogram(dend2,unassignedSamples = c("cluster")))
  expect_silent(plotDendrogram(dend2,leafType="clusters",plotType="colorblock"))
  expect_silent(plotDendrogram(dend2,leafType="samples",plotType="colorblock"))
  
})