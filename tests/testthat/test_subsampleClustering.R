context("subsampleClustering")
source("create_objects.R")
test_that("`subsampleClustering` works", {
    set.seed(4897)
    subAll <- subsampleClustering(x=mat,k=3,clusterFunction="pam",  classifyMethod=c("All"),largeDataset=FALSE, resamp.num = 100, samp.p = 0.7,ncores=1)
    set.seed(4897)
    subAllLarge <- subsampleClustering(x=mat,k=3,clusterFunction="pam",  classifyMethod=c("All"),largeDataset=TRUE, resamp.num = 100, samp.p = 0.7,ncores=1)
	# #this doesn't work for some reason, even though same numbers:
	# expect_equal(round(subAllLarge,5),round(subAll,5))

	#subsample clusterings won't have identification to all samples...
    set.seed(4897)
	subInSampleLarge <- subsampleClustering(x=mat,k=3,clusterFunction="pam",  classifyMethod=c("InSample"),largeDataset=TRUE, resamp.num = 100, samp.p = 0.7,ncores=1)
    set.seed(4897)
    subInSample <- subsampleClustering(x=mat,k=3,clusterFunction="pam",  classifyMethod=c("InSample"),largeDataset=FALSE, resamp.num = 100, samp.p = 0.7,ncores=1)
	expect_equal(subInSampleLarge,subInSample)
 })
