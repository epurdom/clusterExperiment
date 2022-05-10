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

test_that("reduce and filter work with hdf5",{
	expect_silent(filterNames(hdfSCE))

	expect_silent(f1<-filterData(hdfSCE,filterStats="Filter1",cutoff=1))
	expect_silent(f1<-filterData(hdfSCE,filterStats="Filter1",cutoff=1))

	expect_silent(fs<-makeFilterStats(hdfSCE,filterStats="var"))
	expect_silent(fs<-makeFilterStats(hdfSCE,filterStats=c("mean","var")))
	expect_silent(out<-filterData(fs,filterStats=c("mean"),cutoff=1))
	expect_silent(fs<-makeFilterStats(assay(hdfSCE),filterStats="var"))


	expect_silent(defaultNDims(hdfObj,"PCA"))

	#add pca to it
	nDim<-3
	expect_silent(dr3<-makeReducedDims(hdfObj,reducedDims="PCA",maxDims=nDim))
	expect_equal(defaultNDims(dr3,"PCA"),3)

	#test directly on hdf5
	expect_silent(dr3<-makeReducedDims(assay(hdfObj),reducedDims="PCA",maxDims=nDim))


	#test transformation -- need make CE object
    expect_silent(clustNothing1 <- clusterSingle(hdfObj,
		  reduceMethod = "none",
  		  mainClusterArgs=list( clusterArgs=list(k=3),clusterFunction=listBuiltInFunctions()[[1]]),
  	       subsample=FALSE, sequential=FALSE, isCount=FALSE))

	expect_silent(transformation(clustNothing1)<-function(x){exp(x)})
	expect_equal(exp(assay(clustNothing1)),unname(transformData(clustNothing1)))



})

test_that("`getBestFeatures` works with HDF5 assay slot",{
    expect_silent(cl1 <- clusterSingle(hdfObj, 
            subsample=FALSE, sequential=FALSE,
			mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
			isCount=FALSE))
    expect_silent(getBestFeatures(cl1,DEMethod="limma"))
								
	
})