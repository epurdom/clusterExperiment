context("subsampleClustering")
source("create_objects.R")
test_that("`subsampleClustering` works", {
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
	
	#test in passing to subsampleArgs
	rsecOut1<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5, largeDataset=TRUE),random.seed=495)
   	rsecOut2<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
          subsampleArgs=list(resamp.num=5, largeDataset=FALSE),random.seed=495)
	expect_identical(clusterMatrix(rsecOut1),clusterMatrix(rsecOut2))
  
 })
