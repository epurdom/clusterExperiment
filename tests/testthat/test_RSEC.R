context("RSEC")
source("create_objects.R")
test_that("`RSEC` works with matrix, ClusterExperiment, summarizedExperiment",{
	##these examples don't do dendrogram/merge because all -1 after combineMany
	##only tests clusterMany, combineMany parts.
	##so can't do expect_silent, because returns NOTE about that issue.
	rsecOut1<-RSEC(x=mat, isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
   rsecOut2<-RSEC(x=cc, isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  rsecOut3<-RSEC(x=ccSE,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495)
   rsecOut4<-RSEC(x=se,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495)
#test rerunClusterMany argument:
	    rsecOut5<-RSEC(rsecOut2,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",rerunClusterMany=TRUE,subsampleArgs=list(resamp.num=5),random.seed=495)
	   #makes dendrogram so important have here so has to catch defaults of RSEC...
	    rsecOut6<-RSEC(rsecOut2,isCount=FALSE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",rerunClusterMany=FALSE,subsampleArgs=list(resamp.num=5),random.seed=495)
   })

test_that("`RSEC` works through whole series of steps",{
#bigger example where actually goes through all the steps, takes some time:
rsecOut<-RSEC(x=assay(seSimCount), isCount=TRUE,dimReduce="none",
              k0s=4:5,clusterFunction="tight", alphas=0.1,
              betas=0.9,dendroReduce="none",minSizes=1,
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  ##check same as individual steps
 expect_silent(ceOut<-clusterMany(x=assay(seSimCount),ks=4:5,clusterFunction="tight",alphas=0.1,betas=0.9,minSizes=1,
  isCount=TRUE, dimReduce="none", transFun = NULL,
 sequential=TRUE,removeSil=FALSE,subsample=TRUE,silCutoff=0,distFunction=NA,
                 nVarDims=NA,nPCADims=NA,
                 mainClusterArgs=NULL,subsampleArgs=list(resamp.num=5),
                 ncores=1,run=TRUE,seqArgs=list(verbose=FALSE),random.seed=495
 ))
	expect_equal(clusterMatrix(rsecOut,whichClusters="clusterMany"),clusterMatrix(ceOut))
	#gives 'note', can't use expect_silent
 combOut<-combineMany(ceOut, proportion = 0.7,minSize = 5)
 expect_equal(clusterMatrix(rsecOut,whichClusters="combineMany"),clusterMatrix(combOut,whichClusters="combineMany"))
 expect_equal(coClustering(rsecOut),coClustering(combOut))

 expect_silent(dendOut<-makeDendrogram(combOut,dimReduce="none",ndims=NA))
 expect_equal(dendOut@dendro_clusters,rsecOut@dendro_clusters)
 expect_equal(dendOut@dendro_outbranch,rsecOut@dendro_outbranch)

 #now should be the same, check all objects except dendro_samples because very big:
 mergeOut<-mergeClusters(dendOut,mergeMethod = "adjP", cutoff = 0.05,isCount=TRUE)
 expect_equal(dendroClusterIndex(mergeOut),dendroClusterIndex(rsecOut))
 expect_equal(mergeOut@dendro_clusters,rsecOut@dendro_clusters)
 expect_equal(mergeOut@dendro_outbranch,rsecOut@dendro_outbranch)
 expect_equal(coClustering(mergeOut),coClustering(rsecOut))
 expect_equal(clusterMatrix(rsecOut,whichClusters="mergeClusters"), clusterMatrix(mergeOut,whichClusters="mergeClusters"))
 expect_equal(clusterTypes(rsecOut),clusterTypes(mergeOut))
})

test_that("`RSEC` works with no merging",{
  #bigger example where actually goes through all the steps (above skips the merging, in particular, because no dendrogram); takes some time:
  rsecOut<-RSEC(x=assay(seSimCount), isCount=TRUE,dimReduce="none",
                k0s=4:5,clusterFunction="tight", alphas=0.1,
                betas=0.9,dendroReduce="none",minSizes=1,
                subsampleArgs=list(resamp.num=5),random.seed=495,
                mergeMethod="none")
})
