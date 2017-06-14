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
#bigger example where actually goes through all the steps (above skips the merging, in particular, because no dendrogram):
  test<-RSEC(x=seSimCount, isCount=TRUE,dimReduce="none",k0s=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  #test have coClustering object

  #test have dendrogram slots

  ##do individual steps that should be same?
  # ceOut<-clusterMany(x=cc,isCount=FALSE,dimReduce="none",ks=4:5,clusterFunction="tight", alphas=0.1,subsample=TRUE,sequential=TRUE)

})
