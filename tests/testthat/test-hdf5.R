test_that("`ClusterExperiment` constructor works with hdf5",{
    #test creation
    expect_silent(ClusterExperiment(assay(hdfSCE), sample(1:3,size=ncol(hdfSCE),replace=TRUE)))
	
	
})

test_that("`makeDendrogram` works with hdf5",{
    expect_silent(clustNothing <- clusterMany(hdfSCE,
		 ks=c(3,4),clusterFunction="pam",
		 subsample=FALSE, sequential=FALSE,
		 isCount=FALSE,verbose=FALSE))
    expect_silent(clustNothing<-makeConsensus(clustNothing, proportion=1,whichClusters = "clusterMany"))
	expect_silent(makeDendrogram(clustNothing))
	
})
