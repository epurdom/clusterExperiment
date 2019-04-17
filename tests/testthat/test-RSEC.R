context("RSEC")
seedValue<-495 #01875 works for sample.kind="Rejection" (i.e. can merge without error that all -1)
test_that("`RSEC` works with matrix, ClusterExperiment, summarizedExperiment",{
	##these examples don't do dendrogram/merge because all -1 after makeConsensus
	##only tests clusterMany, makeConsensus parts.
	##so can't do expect_silent, because returns NOTE about that issue.
	set.seed(seedValue)
	expect_message(rsecOut1<-RSEC(x=mat, isCount=FALSE,reduceMethod="none",k0s=4:5,
		clusterFunction="tight", alphas=0.1,dendroReduce="none",
        subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE
  	 	),"makeDendrogram encountered following error")
  set.seed(seedValue)
	expect_message(rsecOut2<-RSEC(x=cc, reduceMethod="none",k0s=4:5,
   		clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE
  	 	),"makeDendrogram encountered following error")
  set.seed(seedValue)
	expect_message(rsecOut3<-RSEC(x=ccSE,
	  reduceMethod="none",k0s=4:5,clusterFunction="tight", 
	  alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
	   "makeDendrogram encountered following error")
   set.seed(seedValue)
	 expect_message(rsecOut4<-RSEC(x=se,isCount=FALSE,reduceMethod="none",
   		k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
	   "makeDendrogram encountered following error")
	   #test rerunClusterMany argument:
	set.seed(seedValue)
	expect_message(rsecOut5<-RSEC(rsecOut2,reduceMethod="none",
		k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",rerunClusterMany=TRUE,
		subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
		"makeDendrogram encountered following error")
    #makes dendrogram so important have here so has to catch defaults of RSEC...
	set.seed(seedValue)
	expect_message(rsecOut6<-RSEC(rsecOut2,
			reduceMethod="none",k0s=4:5,clusterFunction="tight", 
			alphas=0.1,dendroReduce="none",rerunClusterMany=FALSE,
			subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
			"makeDendrogram encountered following error")
})

test_that("`RSEC` works through whole series of steps",{
#bigger example where actually goes through all the steps, takes some time:
  set.seed(seedValue)
	expect_message(rsecOut<-RSEC(x=assay(seSimCount), isCount=TRUE,reduceMethod="none",
              k0s=4:5,clusterFunction="tight", alphas=0.1,
              betas=0.9,dendroReduce="none",minSizes=1,
       subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
  "Merging will be done on")
   set.seed(seedValue)
	 expect_silent(ceOut<- clusterMany(x=assay(seSimCount), ks=4:5, clusterFunction="tight", alphas=0.1, betas=0.9, minSizes=1,
  isCount=TRUE, reduceMethod="none", transFun = NULL,
 sequential=TRUE,removeSil=FALSE,subsample=TRUE,silCutoff=0,distFunction=NA,
                 nFilterDims=NA,nReducedDims=NA,
                 mainClusterArgs=NULL,subsampleArgs=list(resamp.num=5),
                 mc.cores=1,run=TRUE,seqArgs=list(verbose=FALSE),mc.set.seed=FALSE
 ))
	expect_equal(clusterMatrix(rsecOut,whichClusters="clusterMany"),clusterMatrix(ceOut))
 expect_message(combOut<-makeConsensus(ceOut, proportion = 0.7,minSize = 5),"no clusters specified to combine")
  expect_equal(clusterMatrix(rsecOut,whichClusters="makeConsensus"),clusterMatrix(combOut,whichClusters="makeConsensus"))
 expect_equal(coClustering(rsecOut),coClustering(combOut))

 expect_silent(dendOut<-makeDendrogram(combOut,reduceMethod="none",nDims=NA))
#they differ in tdata:
 expect_equal(as(dendOut@dendro_clusters,"phylo4"), as(rsecOut@dendro_clusters,"phylo4"))
 tdRsec<-phylobase::tdata(rsecOut@dendro_clusters)
 tdRsec<-tdRsec[,-grep("ClusterIdMerge",colnames(tdRsec))]
 td<-phylobase::tdata(dendOut@dendro_clusters)
 td<-td[,-grep("ClusterIdMerge",colnames(td))]
 expect_equal(td, tdRsec)
 expect_equal(clusterExperiment:::.hasOutBranch(dendOut), clusterExperiment:::.hasOutBranch(rsecOut))
 #now should be the same, check all objects except dendro_samples because very big:
 mergeOut<-mergeClusters(dendOut,mergeMethod = "adjP", DEMethod="edgeR",cutoff = 0.05)
 expect_equal(dendroClusterIndex(mergeOut),dendroClusterIndex(rsecOut))
 expect_equal(mergeOut@dendro_clusters,rsecOut@dendro_clusters)
 expect_equal(clusterExperiment:::.hasOutBranch(mergeOut),clusterExperiment:::.hasOutBranch(rsecOut))
 expect_equal(coClustering(mergeOut),coClustering(rsecOut))
 expect_equal(clusterMatrix(rsecOut,whichClusters="mergeClusters"), clusterMatrix(mergeOut,whichClusters="mergeClusters"))
 expect_equal(clusterTypes(rsecOut),clusterTypes(mergeOut))
})

test_that("`RSEC` works with no merging",{
  #bigger example where actually goes through all the steps (above skips the merging, in particular, because no dendrogram); takes some time:
  set.seed(seedValue)
	rsecOut<-RSEC(x=assay(seSimCount), isCount=TRUE,reduceMethod="none",
                k0s=4:5,clusterFunction="tight", alphas=0.1,
                betas=0.9,dendroReduce="none",minSizes=1,
                subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE,
                mergeMethod="none")
})

test_that("`RSEC` returns clusterMany even when errors later",{
	#error in makeConsensus param
	set.seed(seedValue)
	expect_message(rsecOut1<-RSEC(x=mat, isCount=FALSE,k0s=4:5,
		clusterFunction="tight", alphas=0.1, nReducedDims=3,
        subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE, consensusProportion = -1, consensusMinSize = 5),"Invalid value for the 'proportion' parameter"
  	 	)
	expect_true("clusterMany" %in% clusterTypes(rsecOut1))

	#error in dendro param
	set.seed(seedValue)
	expect_message(rsecOut2<-RSEC(x=mat, isCount=FALSE,k0s=4:5,
		clusterFunction="tight", alphas=0.1, nReducedDims=3,
        subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE, 
		dendroReduce="myfakemethod"
  	 	),"does not contain the given 'reduceMethod' value")
    expect_true(all(c("clusterMany","makeConsensus") %in% clusterTypes(rsecOut2)))
	
	#error in merging -- have to get one where can make dendrogram... takes longer.
	set.seed(seedValue)
	expect_message(rsecOut3<-RSEC(x=assay(seSimCount[sample(size=50,x=1:nrow(seSimCount)),]), isCount=TRUE,reduceMethod="none",
              k0s=4:5,clusterFunction="pam", alphas=0.1,
              betas=0.9,dendroReduce="none",minSizes=1,
       subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE,
		mergeMethod="fakeMerge"
  	 	),"mergeClusters encountered following error")
    expect_true(all(c("clusterMany","makeConsensus") %in% clusterTypes(rsecOut3)))

})

test_that("`RSEC` works with hdf5",{
	#no reduce method, do everything on raw data
	#currently error: Error in tcrossprod(x, y) : 
#  requires numeric/complex matrix/vector arguments

	set.seed(seedValue)
	expect_message(rsecOut1<-RSEC(hdfObj, isCount=FALSE,k0s=4:5,reduceMethod="none",
		clusterFunction="tight", alphas=0.1, 
        subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
		"All samples are unassigned for"
		)

	set.seed(seedValue)
	expect_message(rsecOut2<-RSEC(hdfObj, isCount=FALSE,k0s=4:5,reduceMethod="PCA",
		clusterFunction="tight", alphas=0.1, nReducedDims=3,
        subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
		"Merging will be done on"
		)

	set.seed(seedValue)
	expect_message(rsecOut3<-RSEC(assay(hdfObj), isCount=FALSE,k0s=4:5,reduceMethod="PCA",
		clusterFunction="tight", alphas=0.1, nReducedDims=3,
	    subsampleArgs=list(resamp.num=5),mc.set.seed=FALSE),
		"Merging will be done on"
		)

	expect_equal(clusterMatrix(rsecOut2),clusterMatrix(rsecOut3))
})
