context("makeConsensus")


test_that("`makeConsensus` works with matrix and ClusterExperiment objects", {
    expect_silent(clustNothing <- clusterMany(mat,
        ks=c(3,4),clusterFunction="pam",
        subsample=FALSE, sequential=FALSE,
        isCount=FALSE,verbose=FALSE))
    expect_silent(x1<-makeConsensus(clustNothing,
        proportion=1,whichClusters = "clusterMany"))
    expect_message(x2<-makeConsensus(clustNothing,
        proportion=1),"no clusters specified to combine")
    expect_equal(x1,x2)


    expect_silent(shared1 <- makeConsensus(clusterMatrix(clustNothing),
        proportion=1))
    expect_silent(shared2 <- makeConsensus(clustNothing, 
        "all",proportion=1))
    expect_equal(shared1, primaryCluster(shared2))

    expect_silent(shared3 <- makeConsensus(clustNothing,
         "workflow",proportion=1))
    expect_equal(shared2, shared3)

    expect_silent(shared4 <- makeConsensus(clustNothing,
        1:nClusterings(clustNothing),proportion=1))
    expect_equal(shared3, shared4)

    expect_silent(shared5 <- makeConsensus(clustNothing, 
        "workflow", proportion=.5))

    expect_silent(shared6 <- makeConsensus(clustNothing, 
        "workflow",
        proportion=.5, clusterFunction="tight"))
	expect_true("makeConsensus" %in% clusterTypes(shared6))
	expect_true("makeConsensus" %in% clusterLabels(shared6))
	expect_true(all(primaryCluster(shared6)==
        clusterMatrix(shared6)[,"makeConsensus"]))

  	#----
  	##Check error catching:
  	#----
	expect_silent(ceObj<-clusterSingle(mat, subsample=FALSE,
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3))))
  	expect_error(makeConsensus(ceObj,proportion=1),
                   "no clusters specified")
 	expect_error(makeConsensus(ceObj,
        whichClusters = "all",proportion=-1),
           "Invalid value for the 'proportion' parameter")
           
    #     # FIXME: This should work, once I fix the 'cat' stuff.
    #     # right now just problem of order of when do calculations.
    # expect_error(makeConsensus(ceObj,
    #         whichClusters = "all",proportion=0.7, minSize=-3),
    #        "Invalid value for the 'minSize' parameter")
	expect_error(makeConsensus(ceObj,
        whichClusters = "all",proportion=0.7, propUnassigned=3),
	       "Invalid value for the 'propUnassigned' parameter")
	expect_error(makeConsensus(ceObj,
        whichCluster="clusterSingle"),
  	        'argument "proportion" is missing, with no default')
  	expect_error(makeConsensus(ceObj,
        whichCluster="clusterSingle"),
  			'argument "proportion" is missing, with no default')
    expect_error(makeConsensus(clustNothing, "workflow",
           proportion=.5,clusterFunction="pam"), 
		   "makeConsensus is only implemented for '01' type")

})

test_that("`makeConsensus` works with hdf5",{
    expect_silent(clustNothing <- clusterMany(hdfSCE,
		 ks=c(3,4),clusterFunction="pam",
		 subsample=FALSE, sequential=FALSE,
		 isCount=FALSE,verbose=FALSE))
    expect_silent(x1<-makeConsensus(clustNothing,proportion=1,whichClusters = "clusterMany"))
	
})

test_that("`makeConsensus` works when multiple runs of workflow", {
  expect_silent(clustNothing <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE))

  expect_silent(shared1 <- makeConsensus(clustNothing, "all",proportion=1))
  expect_silent(shared2<-makeConsensus(shared1,"all",proportion=1))
  expect_true("makeConsensus.1" %in% clusterTypes(shared2))

  expect_silent(clustNothing2 <- clusterMany(shared2, ks=c(5,6), clusterFunction="pam",
                               subsample=FALSE, sequential=FALSE, verbose=FALSE))
  expect_true("makeConsensus.1" %in% clusterTypes(clustNothing2))
  expect_true("clusterMany.2" %in% clusterTypes(clustNothing2))
  expect_true("makeConsensus.2" %in% clusterTypes(clustNothing2))

  expect_silent(shared3 <- makeConsensus(clustNothing2, "all",proportion=1))
  expect_silent(shared4 <- makeConsensus(clusterMatrix(clustNothing2),proportion=1))
  expect_equal(shared4, primaryCluster(shared3))

  expect_silent(shared5 <- makeConsensus(clustNothing2, "workflow",proportion=1))
  expect_silent(shared6 <- makeConsensus(clusterMatrix(clustNothing2)[,1:2],proportion=1))
  expect_equal(shared6, primaryCluster(shared5))


  expect_silent(clustNothing3 <- addClusterings(clustNothing2, primaryCluster(shared5)))
  expect_silent(shared7 <- makeConsensus(clustNothing3, "all",proportion=1))
  expect_silent(shared8 <- makeConsensus(clustNothing3, "workflow",proportion=1))
})


test_that("`makeConsensus` more esoteric options", {
    expect_silent(clustNothing <- clusterMany(mat,
       ks=c(3,4),clusterFunction="pam",
       subsample=FALSE, sequential=FALSE,
       isCount=FALSE,verbose=FALSE))

    #----
    #check passing to clusterArgs
    #----
    expect_silent(shared1 <- makeConsensus(clustNothing, "all",
        proportion=0.7,clusterArgs=list(removeDup=FALSE)))
    #silently ignores passing alpha
    expect_silent(shared2 <- makeConsensus(clustNothing, "all",
        proportion=0.7,clusterArgs=list(removeDup=FALSE,alpha=0.5)))
    expect_equal(primaryCluster(shared1),primaryCluster(shared2))
    expect_silent(makeConsensus(clustNothing, "all",
        proportion=0.7,clusterArgs=list(removeDup=FALSE,
            evalClusterMethod=c("maximum"))))
    
    #----
    # Check whenUnassign
    #----
    expect_silent(makeConsensus(clustNothing, "all",
        proportion=0.7,whenUnassign="before"))
    expect_silent(makeConsensus(clustNothing, "all",
        proportion=0.7,whenUnassign="after"))

})

test_that("`RSEC` pass options to `makeConsensus`", {
	## these examples don't do dendrogram/merge because all -1 after makeConsensus
	## only tests clusterMany, makeConsensus parts.
	## so can't do expect_silent, because returns NOTE about that issue.
    seedValue<-495 
	expect_message(rsecOut1<-RSEC(x=mat,
        isCount=FALSE,reduceMethod="none",k0s=4:5,
		clusterFunction="tight", alphas=0.1,dendroReduce="none",
        subsampleArgs=list(resamp.num=5), consensusProportion=0.7,
        consensusMinSize=3, sequential=FALSE,
        consensusArgs=list(clusterFunction="tight"),
        random.seed=seedValue
  	 	),"makeDendrogram encountered following error")
    test<-makeConsensus(clusterMatrix(rsecOut1,
        whichClusters="clusterMany"),
        proportion=0.7,minSize=3,
        clusterFunction="tight")
    expect_equal(test,
        clusterMatrix(rsecOut1,"makeConsensus")[,1])

	expect_message(rsecOut1<-RSEC(x=mat,
        isCount=FALSE,reduceMethod="none",k0s=4:5,
		clusterFunction="tight", alphas=0.1,dendroReduce="none",
        subsampleArgs=list(resamp.num=5),
        consensusProportion=0.7,sequential=FALSE,
        consensusMinSize=3,
        consensusArgs=list(
            clusterFunction="hierarchical01",
            clusterArgs=list(removeDup=FALSE,
                evalClusterMethod=c("maximum") )),
        random.seed=seedValue
  	 	),"makeDendrogram encountered following error")
    test<-makeConsensus(clusterMatrix(rsecOut1,
        whichClusters="clusterMany"),
        proportion=0.7,minSize=3,
        clusterArgs=list(removeDup=FALSE,
            evalClusterMethod=c("maximum")),
        clusterFunction="hierarchical01")
    expect_equal(test,
        clusterMatrix(rsecOut1,"makeConsensus")[,1])

})

test_that("`makeConsensus` preserves the colData and rowData of SE", {
    expect_silent(cl <- clusterMany(se, verbose=FALSE,
        clusterFunction="pam", ks=2:4, findBestK=c(FALSE),
        removeSil=TRUE, subsample=FALSE))

    expect_silent(cl <- makeConsensus(cl, 
        whichClusters="workflow", proportion=0.5))

    expect_equal(colData(cl),colData(se))
    expect_equal(rownames(cl),rownames(se))
    expect_equal(colnames(cl),colnames(se))
    expect_equal(metadata(cl),metadata(se))
    expect_equal(rowData(cl),rowData(se))

})
