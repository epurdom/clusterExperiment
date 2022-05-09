test_that("`makeDendrogram` works with RSECClass objects", {
    expect_silent(cl1 <- clusterSingle(smSimData, subsample=FALSE, sequential=FALSE, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),isCount=FALSE))
    expect_silent(clustWithDendro <- makeDendrogram(cl1))
	#--------
    #Check remove clusters when have dendrogram
	#--------
	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adj", plotInfo="none",DEMethod="limma"),"Merging will be done on ' clusterSingle ', with clustering index 1")
    expect_silent(removeClusterings(clustMerged,whichClusters="mergeClusters")) #remove merged, keep one with dendrogram
    expect_silent(removeClusterings(clustMerged,whichClusters=2)) #remove one with dendrogram

    bigCE<-ceSimCount
    expect_silent(bigCE<-makeDendrogram(bigCE,whichCluster="cluster1"))
    x1<-bigCE
    #--- check clusterMany updates dendrogram correctly
    expect_silent(bigCE<-clusterMany(bigCE,
        k=2:8,clusterFunction="hierarchicalK",
        verbose=FALSE,
        makeMissingDiss=TRUE))
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],
        clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    #takes a long time!
    #expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    expect_error(makeDendrogram(bigCE,whichCluster="workflow"),
        "Invalid value for 'whichCluster'")
 
    #--- check makeConsensus updates dendrogram correctly
    expect_message(bigCE<-makeConsensus(bigCE,proportion=0.3),
        "no clusters specified to combine, using results from clusterMany")
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],
        clusterLabels(x1)[x1@dendro_index])
    expect_equal(bigCE@dendro_clusters,x1@dendro_clusters) 
    expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    expect_silent(makeDendrogram(bigCE,whichCluster="makeConsensus") )
    
    
    #--- check mergeClusters updates dendrogram correctly
    expect_message(bigCE<-mergeClusters(bigCE,
        mergeMethod="adjP",cutoff=0.2, DEMethod="limma"),
        "Merging will be done on ")
    expect_equal(clusterLabels(bigCE)[bigCE@dendro_index],
        clusterLabels(x1)[x1@dendro_index])
    #Note, x1 doesn't give any merged clusters in dendrogram because before mergeClusters step...  so going to remove that element from bigCE
	tdf<-phylobase::tdata(bigCE@dendro_clusters)
	tdf$ClusterIdMerge<-NA
	phylobase::tdata(bigCE@dendro_clusters)<-tdf
	expect_equal(bigCE@dendro_clusters,x1@dendro_clusters)  
    expect_equal(bigCE@dendro_samples,x1@dendro_samples) 
    
    expect_error(getBestFeatures(bigCE,contrastType="Dendro"),
        "only single cluster in clustering -- cannot run getBestFeatures")
	expect_silent(primaryClusterIndex(bigCE)<-3)
	expect_error( getBestFeatures(bigCE,contrastType="Dendro"),"does not match either the cluster on which the dendrogram was made or the merge cluster from this dendrogram")
}