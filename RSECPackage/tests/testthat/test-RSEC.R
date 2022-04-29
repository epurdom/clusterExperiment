context("RSEC")
seedValue<-495 #01875 works for sample.kind="Rejection"
test_that("`RSEC` works with matrix, ClusterExperiment, summarizedExperiment",{
	## these examples don't do dendrogram/merge because all -1 after makeConsensus
	## only tests clusterMany, makeConsensus parts.
	## so can't do expect_silent, because returns NOTE about that issue.
	expect_message(rsecOut1<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
		clusterFunction="tight", alphas=0.1,dendroReduce="none",
        subsampleArgs=list(resamp.num=5),random.seed=seedValue
  	 	),"makeDendrogram encountered following error")
   expect_message(rsecOut2<-RSEC(x=cc, reduceMethod="none",k0s=4:5,
   		clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=seedValue
  	 	),"makeDendrogram encountered following error")
  expect_message(rsecOut3<-RSEC(x=ccSE,
	  reduceMethod="none",k0s=4:5,clusterFunction="tight", 
	  alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=seedValue),
	   "makeDendrogram encountered following error")
   expect_message(rsecOut4<-RSEC(x=se,isCount=FALSE,reduceMethod="none",
   		k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=seedValue),
	   "makeDendrogram encountered following error")
	   #test rerunClusterMany argument:
	expect_message(rsecOut5<-RSEC(rsecOut2,reduceMethod="none",
		k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",rerunClusterMany=TRUE,
		subsampleArgs=list(resamp.num=5),random.seed=seedValue),
		"makeDendrogram encountered following error")
    #makes dendrogram so important have here so has to catch defaults of RSEC...
	expect_message(rsecOut6<-RSEC(rsecOut2,
			reduceMethod="none",k0s=4:5,clusterFunction="tight", 
			alphas=0.1,dendroReduce="none",rerunClusterMany=FALSE,
			subsampleArgs=list(resamp.num=5),random.seed=seedValue),
			"makeDendrogram encountered following error")
})



test_that("`RSEC` works with no merging",{
  #do the same, only don't do merging:
    expect_message(rsecOut<-RSEC(x=assay(seSimCount), 
        isCount=TRUE,reduceMethod="none",
        k0s=4:5,clusterFunction="hierarchical01", alphas=0.1,
        betas=0.9,dendroReduce="none",minSizes=1,
        seqArgs=list(top.can=0),
        subsampleArgs=list(resamp.num=5),random.seed=seedValue,
        mergeMethod="none"),
        "clusters will not be merged because argument")
})

test_that("`RSEC` returns clusterMany even when errors later",{
	#error in makeConsensus param
	expect_message(rsecOut1<-
        RSEC(x=mat, isCount=FALSE,k0s=4:5,
		clusterFunction="tight", alphas=0.1, nReducedDims=3,
        subsampleArgs=list(resamp.num=5),random.seed=seedValue, 
        consensusProportion = -1, consensusMinSize = 5),
        "Invalid value for the 'proportion' parameter"
  	 	)
	expect_true("clusterMany" %in% clusterTypes(rsecOut1))

	#error in dendro param
	expect_message(rsecOut2<-RSEC(x=mat, isCount=FALSE,k0s=4:5,
		clusterFunction="tight", alphas=0.1, nReducedDims=3,
        subsampleArgs=list(resamp.num=5),random.seed=seedValue, 
		dendroReduce="myfakemethod"
  	 	),"does not contain the given 'reduceMethod' value")
    expect_true(all(c("clusterMany","makeConsensus") 
        %in% clusterTypes(rsecOut2)))
	
	# error in merging -- give wrong method
    # have to get dataset where can actually make dendrogram... 
    # takes longer.
    # > sample(size=50,x=1:nrow(seSimCount))
#      [1]  67  47  34 124  23  40 104  93  18  26 123  77  64  35  90  58 103 153  62  66 138  57 118  50  60 145  52 134  45  78  12
#     [32]  38 121  63 128 111 117  56  43  51 149   9  73 150   6  39 125 133 127 152
    set.seed(4289)
	expect_message(rsecOut3<-RSEC(x=
        assay(seSimCount[sample(size=50,x=1:nrow(seSimCount)),]), 
        isCount=TRUE,reduceMethod="none",
        k0s=4:5,clusterFunction="pam", alphas=0.1,
        betas=0.9,dendroReduce="none",minSizes=1,
        subsampleArgs=list(resamp.num=5,clusterFunction="kmeans"),
        random.seed=seedValue,
		mergeMethod="fakeMerge"
  	 	),"mergeClusters encountered following error")

    expect_true(all(c("clusterMany","makeConsensus") %in% clusterTypes(rsecOut3)))

})


test_that("`RSEC` passing args to subsample",{
    #test in passing to subsampleArgs
  	expect_message(rsecOut1<-RSEC(x=mat, 
  				   isCount=FALSE,reduceMethod="none",k0s=4:5,
  	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
  	               subsampleArgs=list(resamp.num=5),random.seed=495),
  				   "makeDendrogram encountered following error")
  	expect_message(rsecOut2<-RSEC(x=mat, 
  				   isCount=FALSE,reduceMethod="none",k0s=4:5,
  	               clusterFunction="tight", alphas=0.1,dendroReduce="none",
  	               subsampleArgs=list(resamp.num=5),random.seed=495),
  				   "makeDendrogram encountered following error")
  	expect_identical(clusterMatrix(rsecOut1),clusterMatrix(rsecOut2))
    
})
