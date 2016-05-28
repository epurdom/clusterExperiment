context("RSEC")
source("create_objects.R")
test_that("`RSEC` works with matrix, clusterExperiment, summarizedExperiment",{
  ##these examples don't do dendrogram/merge because all -1 after combineMany
  RSEC(x=mat, isCount=FALSE,dimReduce="none",ks=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  RSEC(x=cc, isCount=FALSE,dimReduce="none",ks=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  RSEC(x=ccSE,isCount=FALSE,dimReduce="none",ks=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495)
  RSEC(x=se,isCount=FALSE,dimReduce="none",ks=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495)
#bigger example where actually goes through all the steps:
  RSEC(x=seSimCount, isCount=TRUE,dimReduce="none",ks=4:5,clusterFunction="tight", alphas=0.1,dendroReduce="none",
       subsampleArgs=list(resamp.num=5),random.seed=495
  )
  
})
