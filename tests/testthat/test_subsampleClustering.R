context("subsampleClustering")

test_that("subsampling", {

  set.seed(4897)
  subAll <- subsampleClustering(x=mat,clusterFunction="pam",
                                clusterArgs=list(k=3), classifyMethod="All",
                                resamp.num = 100,
                                samp.p = 0.7,ncores=1)


  ## Windows does not support mclapply
  if(.Platform$OS.type == "unix"){

    set.seed(4897)
    subAllParal <- subsampleClustering(x=mat,clusterFunction="pam",
                                       clusterArgs=list(k=3), classifyMethod="All",
                                       resamp.num = 100,
                                       samp.p = 0.7,ncores=2)

    expect_identical(subAllParal,subAll)

  }


	#subsample clusterings won't have identification to all samples...
  set.seed(4897)
  subInSample <- subsampleClustering(x=mat,clusterArgs=list(k=3),
                                     clusterFunction="pam",
                                     classifyMethod=c("InSample"),
                                     resamp.num = 100,
                                     samp.p = 0.7,ncores=1)

  #test in passing to subsampleArgs
	rsecOut1<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
	               subsampleArgs=list(resamp.num=5),random.seed=495)
	rsecOut2<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
	               subsampleArgs=list(resamp.num=5),random.seed=495)
	expect_identical(clusterMatrix(rsecOut1),clusterMatrix(rsecOut2))

})

