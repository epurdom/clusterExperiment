context("subsampleClustering")

test_that("Large memory option", {

  ## largeDataset = FALSE
  set.seed(4897)
  subAll <- subsampleClustering(x=mat,clusterFunction="pam",
                                clusterArgs=list(k=3), classifyMethod="All",
                                largeDataset=FALSE, resamp.num = 100,
                                samp.p = 0.7,ncores=1)

  ## largeDataset = FALSE; whichImplementation = "R"
  set.seed(4897)
  subAllLarge <- subsampleClustering(x=mat,clusterFunction="pam",
                                     clusterArgs=list(k=3), classifyMethod="All",
                                     largeDataset=TRUE, resamp.num = 100,
                                     samp.p = 0.7,ncores=1)

  ## largeDataset = FALSE; whichImplementation = "Csimple"
  set.seed(4897)
  subAllLarge2 <- subsampleClustering(x=mat,clusterFunction="pam",
                                      clusterArgs=list(k=3), classifyMethod="All",
                                      largeDataset=TRUE, resamp.num = 100,
                                      samp.p = 0.7,ncores=1,
                                      whichImplementation = "Csimple")

  ## largeDataset = FALSE; whichImplementation = "Cmemory"
  set.seed(4897)
  subAllLarge3<- subsampleClustering(x=mat,clusterFunction="pam",
                                     clusterArgs=list(k=3), classifyMethod="All",
                                     largeDataset=TRUE, resamp.num = 100,
                                     samp.p = 0.7,ncores=1,
                                     whichImplementation = "Cmemory")

  expect_identical(subAllLarge,subAll)
  expect_identical(subAllLarge2,subAll)
  expect_identical(subAllLarge3,subAll)

  ## Windows does not support mclapply
  if(.Platform$OS.type == "unix"){

    set.seed(4897)
    subAllParal <- subsampleClustering(x=mat,clusterFunction="pam",
                                       clusterArgs=list(k=3), classifyMethod="All",
                                       largeDataset=FALSE, resamp.num = 100,
                                       samp.p = 0.7,ncores=2)

    set.seed(4897)
    subAllLargeParal <- subsampleClustering(x=mat,clusterFunction="pam",
                                            clusterArgs=list(k=3), classifyMethod="All",
                                            largeDataset=TRUE, resamp.num = 100,
                                            samp.p = 0.7,ncores=2)

    set.seed(4897)
    subAllLarge2Paral <- subsampleClustering(x=mat,clusterFunction="pam",
                                             clusterArgs=list(k=3), classifyMethod="All",
                                             largeDataset=TRUE, resamp.num = 100,
                                             samp.p = 0.7,ncores=2,
                                             whichImplementation = "Csimple")

    set.seed(4897)
    subAllLarge3Paral <- subsampleClustering(x=mat,clusterFunction="pam",
                                             clusterArgs=list(k=3), classifyMethod="All",
                                             largeDataset=TRUE, resamp.num = 100,
                                             samp.p = 0.7,ncores=2,
                                             whichImplementation = "Cmemory")

    expect_identical(subAllLargeParal,subAll)
    expect_identical(subAllLarge2Paral,subAll)
    expect_identical(subAllLarge3Paral,subAll)

  }


	#subsample clusterings won't have identification to all samples...
    set.seed(4897)
	subInSampleLarge <- subsampleClustering(x=mat,clusterArgs=list(k=3),
	                                        clusterFunction="pam",
	                                        classifyMethod=c("InSample"),
	                                        largeDataset=TRUE,
	                                        resamp.num = 100,
	                                        samp.p = 0.7,ncores=1)
    set.seed(4897)
    subInSample <- subsampleClustering(x=mat,clusterArgs=list(k=3),
                                       clusterFunction="pam",
                                       classifyMethod=c("InSample"),
                                       largeDataset=FALSE,
                                       resamp.num = 100,
                                       samp.p = 0.7,ncores=1)
	expect_identical(subInSampleLarge,subInSample)

	#test in passing to subsampleArgs
	rsecOut1<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
	               subsampleArgs=list(resamp.num=5, largeDataset=TRUE),random.seed=495)
	rsecOut2<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
	               subsampleArgs=list(resamp.num=5, largeDataset=FALSE),random.seed=495)
	expect_identical(clusterMatrix(rsecOut1),clusterMatrix(rsecOut2))

	rsecOut3<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
	               subsampleArgs=list(resamp.num=5, largeDataset=TRUE,
	                                  whichImplementation = "Csimple"),random.seed=495)
	expect_identical(clusterMatrix(rsecOut1),clusterMatrix(rsecOut3))

	rsecOut4<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
	               subsampleArgs=list(resamp.num=5, largeDataset=TRUE,
	                                  whichImplementation = "Cmemory"),random.seed=495)
	expect_identical(clusterMatrix(rsecOut1),clusterMatrix(rsecOut4))

 })

