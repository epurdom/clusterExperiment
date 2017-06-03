context("RSEC")
source("create_objects.R")
test_that("`RSEC` works with matrix, clusterExperiment, summarizedExperiment",{
  ##these examples don't do dendrogram/merge because all -1 after combineMany
  ##only tests clusterMany, combineMany parts.
  RSEC(x=mat, isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  rsecOut<-RSEC(x=cc, isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  RSEC(x=ccSE,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495)
  RSEC(x=se,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495)
#test rerunClusterMany argument:
	   RSEC(rsecOut,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",rerunClusterMany=TRUE,subsampleArgs=list(resamp.num=5),random.seed=495)
	   RSEC(rsecOut,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",rerunClusterMany=FALSE,subsampleArgs=list(resamp.num=5),random.seed=495)
   })
  
test_that("`RSEC` works through whole series of steps",{
#bigger example where actually goes through all the steps (above skips the merging, in particular, because no dendrogram); takes some time:
rsecOut<-RSEC(x=assay(seSimCount), isCount=TRUE,dimReduce="none",
              k0s=4:5,clusterFunction="tight", alphas=0.1,
              betas=0.9,dendroReduce="none",minSizes=1,
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  ##check same as individual steps
 ceOut<-clusterMany(x=assay(seSimCount),ks=4:5,clusterFunction="tight",alphas=0.1,betas=0.9,minSizes=1,
  isCount=TRUE, dimReduce="none", transFun = NULL,
 sequential=TRUE,removeSil=FALSE,subsample=TRUE,silCutoff=0,distFunction=NA,
                 nVarDims=NA,nPCADims=NA,
                 clusterDArgs=NULL,subsampleArgs=list(resamp.num=5),
                 ncores=1,run=TRUE,seqArgs=list(verbose=FALSE),random.seed=495
 )
	expect_equal(clusterMatrix(rsecOut,whichClusters="clusterMany"),clusterMatrix(ceOut))

 combOut<-combineMany(ceOut, proportion = 0.7,minSize = 5)
 expect_equal(clusterMatrix(rsecOut,whichClusters="combineMany"),clusterMatrix(combOut,whichClusters="combineMany"))
 expect_equal(coClustering(rsecOut),coClustering(combOut))
 
 dendOut<-makeDendrogram(combOut,dimReduce="none",ndims=NA)
 expect_equal(dendOut@dendro_clusters,rsecOut@dendro_clusters)
 expect_equal(dendOut@dendro_outbranch,rsecOut@dendro_outbranch)
 
 #now should be the same:
 mergeOut<-mergeClusters(dendOut,mergeMethod = "adjP", cutoff = 0.05,isCount=TRUE)
 expect_equal(dendroClusterIndex(mergeOut),dendroClusterIndex(rsecOut))
 expect_equal(clusterMatrix(rsecOut,whichClusters="mergeClusters"),clusterMatrix(mergeOut,whichClusters="mergeClusters"))
})

#code in RSEC:
# ce<-clusterMany(x,ks=k0s,clusterFunction=clusterFunction,alphas=alphas,betas=betas,minSizes=minSizes,
#                 sequential=TRUE,removeSil=FALSE,subsample=TRUE,silCutoff=0,distFunction=NA,
#                 isCount=isCount,transFun=transFun,
#                 dimReduce=dimReduce,nVarDims=nVarDims,nPCADims=nPCADims,
#                 clusterDArgs=clusterDArgs,subsampleArgs=subsampleArgs,
#                 seqArgs=seqArgs,ncores=ncores,random.seed=random.seed,run=run)
#  RSEC(x, isCount = FALSE, transFun = NULL,
#    dimReduce = "PCA", nVarDims = NA, nPCADims = c(50), k0s = 4:15,
#    clusterFunction = c("tight", "hierarchical01"), alphas = c(0.1, 0.2, 0.3),
#    betas = 0.9, minSizes = 1, combineProportion = 0.7,
#    combineMinSize = 5, dendroReduce = "mad", dendroNDims = 1000,
#    mergeMethod = "adjP", mergeCutoff = 0.05, verbose = FALSE,
#    clusterDArgs = NULL, subsampleArgs = NULL, seqArgs = NULL, ncores = 1,
#    random.seed = NULL, run = TRUE)

# 	if("combinedProportion" %in% names(passedArgs)) args1<-c(args1,"proportion"=passedArgs$combineProportion)
	# if("combineMinSize" %in% names(passedArgs)) args1<-c(args1,"minSize"=passedArgs$combineMinSize)
	#   ce<-do.call("combineMany",c(list(x=ce,whichClusters="clusterMany"),args1))


			     #test have coClustering object
	
			     #test have dendrogram slots


#  RSEC(x, isCount = FALSE, transFun = NULL,
#    dimReduce = "PCA", nVarDims = NA, nPCADims = c(50), k0s = 4:15,
#    clusterFunction = c("tight", "hierarchical01"), alphas = c(0.1, 0.2, 0.3),
#    betas = 0.9, minSizes = 1, combineProportion = 0.7,
#    combineMinSize = 5, dendroReduce = "mad", dendroNDims = 1000,
#    mergeMethod = "adjP", mergeCutoff = 0.05, verbose = FALSE,
#    clusterDArgs = NULL, subsampleArgs = NULL, seqArgs = NULL, ncores = 1,
#    random.seed = NULL, run = TRUE)
#
