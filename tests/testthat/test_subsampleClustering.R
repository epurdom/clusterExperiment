context("subsampleClustering")
source("create_objects.R")
test_that("`subsampleClustering` works", {
	' coOccur <- subsampleClustering(clusterFunction="kmeans", x=simData, 
	#' clusterArgs=list(k=3,nstart=10), resamp.n=100, samp.p=0.7)
	
    set.seed(4897)
    subAll <- subsampleClustering(x=mat,clusterFunction="pam", 
		clusterArgs=list(k=3), classifyMethod="All",
		largeDataset=FALSE, resamp.num = 100, samp.p = 0.7,ncores=1)
    set.seed(4897)
    subAllLarge <- subsampleClustering(x=mat,clusterFunction="pam", 
		clusterArgs=list(k=3), classifyMethod="All",
		largeDataset=TRUE, resamp.num = 100, samp.p = 0.7,ncores=1)
	# #this doesn't work for some reason, even though same numbers:
    expect_identical(subAllLarge,subAll)

	#subsample clusterings won't have identification to all samples...
    set.seed(4897)
	subInSampleLarge <- subsampleClustering(x=mat,clusterArgs=list(k=3), clusterFunction="pam",  classifyMethod=c("InSample"),largeDataset=TRUE, resamp.num = 100, samp.p = 0.7,ncores=1)
    set.seed(4897)
    subInSample <- subsampleClustering(x=mat,clusterArgs=list(k=3), clusterFunction="pam",  classifyMethod=c("InSample"),largeDataset=FALSE, resamp.num = 100, samp.p = 0.7,ncores=1)
	expect_identical(subInSampleLarge,subInSample)
 })
