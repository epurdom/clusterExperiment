test_that("`ClusterExperiment` constructor works with hdf5",{
    #test creation
    expect_silent(ClusterExperiment(hdfSCE, sample(1:3,size=ncol(hdfSCE),replace=TRUE)))
	
	
})

test_that("`makeDendrogram` works with hdf5",{
  expect_silent(clustNothing<-ClusterExperiment(hdfSCE, sample(1:3,size=ncol(hdfSCE),replace=TRUE)))
    
	expect_silent(makeDendrogram(clustNothing))
	
})

test_that("ClusterExperiment works with hdf5",{
	# Check that make CE keeps reduced dims
  expect_silent(clustNothing1<-ClusterExperiment(hdfSCE, sample(1:3,size=ncol(hdfSCE),replace=TRUE)))
	expect_equal(reducedDimNames(hdfSCE), reducedDimNames(clustNothing1))
	 
	#test transformation (on CE object)
	expect_silent(transformation(clustNothing1)<-function(x){exp(x)})
	expect_equal(exp(assay(clustNothing1)),unname(transformData(clustNothing1)))
	
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

	#test directly on hdf5 (not a Summarized experiment)
	expect_silent(dr3<-makeReducedDims(assay(hdfObj),reducedDims="PCA",maxDims=nDim))


})

test_that("`getBestFeatures` works with HDF5 assay slot",{
  	expect_silent(clustNothing<-ClusterExperiment(hdfSCE, sample(1:3,size=ncol(hdfSCE),replace=TRUE)))
  
    expect_silent(getBestFeatures(clustNothing,DEMethod="limma"))
								
	
})

test_that("plotting works with hdf5 assays objects",{
	##plotClusters
	expect_silent(cl1<-ClusterExperiment(hdfSCE, sample(1:3,size=ncol(hdfSCE),replace=TRUE)))
	
			#     expect_silent(cl1 <- clusterSingle(hdfSCE, reduceMethod="PCA",
			#             subsample=FALSE, sequential=FALSE,
			# mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
			# isCount=FALSE))
	expect_silent(plotClusters(cl1))
	
	##plotBarplot
	expect_silent(plotBarplot(cl1))
	
	##plotReducedDims
	expect_silent(plotReducedDims(cl1,legend="bottomright"))

	##plotFeatureBoxplot
	expect_silent(plotFeatureBoxplot(object=cl1,feature=1))

  expect_silent(plotHeatmap(hdfObj))
	

})

