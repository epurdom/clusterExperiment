## Construct ClusterExperiment with multiple named assays
suppressMessages(multi_se <- SummarizedExperiment(assays = list(counts = simCount,
                                               logcounts = log1p(simCount))))
suppressMessages(multi_cc <- ClusterExperiment(multi_se, trueCluster))
suppressMessages(ccTrue <- ClusterExperiment(simCount, trueCluster))
suppressMessages(ccTrue2 <- ClusterExperiment(simCount, trueCluster, transformation = log1p))
seedValue<-01875 #works for sample.kind="Rejection"

test_that("RSEC works wih non default assays", {
  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),   
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "counts"),
                 "Merging will be done on")
  
  expect_message(out2<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")
  
  expect_message(out3<-RSEC(x=ccTrue, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue),
                 "Merging will be done on")
  
  expect_message(out4<-RSEC(x=ccTrue2, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue),
                 "Merging will be done on")
  
  expect_false(all(primaryCluster(out1) == primaryCluster(out2)))
  
  expect_equivalent(out1, out3)
  expect_equivalent(out2, out4)
})

test_that("RSEC works independent of assay order", {
  local_edition(3) # allows to ignore function environment; for some reason this wasn't problem before, but now it is
  #create two with different order of the assays
  multi_se <- SummarizedExperiment(assays = list(counts = simCount,
                                                 logcounts = log1p(simCount)))
  multi_se2 <- SummarizedExperiment(assays = list(logcounts = log1p(simCount),
                                                  counts = simCount))
  multi_cc <- ClusterExperiment(multi_se, trueCluster)
  multi_cc2 <- ClusterExperiment(multi_se2, trueCluster)
  
  #use character, logcounts on both 
  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")
  
  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")
  
  # had to change the tests, because with new version of bioconductor, 
  # order of assays being different is caught, and for some reason wasn't caught before. 
  # Very strange passed previously.
  # Want to test everything same but assays
  out1@assays<-out2@assays
  expect_equal(out1, out2,ignore_function_env=TRUE)
  
  #use numeric
  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = 2),
                 "Merging will be done on")
  
  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="none",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = 1),
                 "Merging will be done on")
  # Want to test everything same but assays
  out1@assays<-out2@assays
  expect_equal(out1, out2,ignore_function_env=TRUE)
  
  #use character, counts on both with PCA reduce
  # In 2.7.2-9001 had to change this to 60 dims -- no longer ran, probably because changed the clustering slightly?
  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="PCA", nReducedDims = 60,
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5, clusterFunction="pam"),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "counts"),
                 "Merging will be done on")
  
  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="PCA", nReducedDims = 60,
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5, clusterFunction="pam"),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "counts"),
                 "Merging will be done on")
  
  # Want to test everything same but assays
  out1@assays<-out2@assays
  expect_equal(out1, out2,ignore_function_env=TRUE)
  
  #use character, counts on both with var reduce
  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="var",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")
  
  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="var",
                            k0s=4:5, clusterFunction="hierarchical01", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            seqArgs=list(top.can=0),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")
  
  # Want to test everything same but assays
  out1@assays<-out2@assays
  expect_equal(out1, out2,ignore_function_env=TRUE)
  
})

test_that("`RSEC` works through whole series of steps",{
  #bigger example where actually goes through all the steps, takes some time:
  expect_message(rsecOut<-RSEC(
    x=assay(seSimCount), isCount=TRUE,reduceMethod="none",
    k0s=4:5,clusterFunction="hierarchical01", alphas=0.1, 
    betas=0.9,minSizes=1,
    subsampleArgs=list(resamp.num=5), seqArgs=list(top.can=0),
    random.seed=seedValue, makeMissingDiss=TRUE,
    consensusProportion=0.7, consensusMinSize=5,
    dendroReduce="none", stopOnErrors = FALSE,
    mergeMethod = "adjP", 
    mergeDEMethod="edgeR",mergeCutoff = 0.05),
    "Merging will be done on")
  expect_silent(ceOut<- clusterMany(x=assay(seSimCount), isCount=TRUE,
                                    reduceMethod="none", 
                                    ks=4:5, clusterFunction="hierarchical01", alphas=0.1, 
                                    betas=0.9, minSizes=1,
                                    transFun = NULL,
                                    sequential=TRUE,removeSil=FALSE,subsample=TRUE,
                                    silCutoff=0,distFunction=NA,
                                    nFilterDims=NA,nReducedDims=NA,
                                    mainClusterArgs=NULL,subsampleArgs=list(resamp.num=5),
                                    ncores=1,run=TRUE, verbose=FALSE,stopOnErrors = TRUE,
                                    seqArgs=list(verbose=FALSE,top.can=0),random.seed=seedValue
  ))
  expect_equal(clusterMatrix(rsecOut,
                             whichClusters="clusterMany"),clusterMatrix(ceOut))
  expect_message(combOut<-makeConsensus(ceOut, 
                                        proportion = 0.7,minSize = 5),
                 "no clusters specified to combine")
  expect_equal(clusterMatrix(rsecOut,
                             whichClusters="makeConsensus"),
               clusterMatrix(combOut,
                             whichClusters="makeConsensus"))
  expect_equal(clusterMatrix(rsecOut,which=coClustering(rsecOut)),
               clusterMatrix(combOut,which=coClustering(combOut)))
  
  expect_silent(dendOut<-makeDendrogram(combOut,
                                        reduceMethod="none",nDims=NA))
  #they differ in tdata:
  expect_equal(as(dendOut@dendro_clusters,"phylo4"), 
               as(rsecOut@dendro_clusters,"phylo4"))
  tdRsec<-phylobase::tdata(rsecOut@dendro_clusters)
  tdRsec<-tdRsec[,-grep("ClusterIdMerge",colnames(tdRsec))]
  td<-phylobase::tdata(dendOut@dendro_clusters)
  td<-td[,-grep("ClusterIdMerge",colnames(td))]
  expect_equal(td, tdRsec)
  expect_equal(clusterExperiment:::.hasOutBranch(dendOut),
               clusterExperiment:::.hasOutBranch(rsecOut))
  #now should be the same, check all objects except dendro_samples because very big:
  expect_message(mergeOut<-mergeClusters(dendOut,mergeMethod = "adjP", 
                                         DEMethod="edgeR",cutoff = 0.05),
                 "Merging will be done on")
  expect_equal(dendroClusterIndex(mergeOut),dendroClusterIndex(rsecOut))
  expect_equal(mergeOut@dendro_clusters,rsecOut@dendro_clusters)
  expect_equal(clusterExperiment:::.hasOutBranch(mergeOut),
               clusterExperiment:::.hasOutBranch(rsecOut))
  expect_equal(coClustering(mergeOut),coClustering(rsecOut))
  expect_equal(clusterMatrix(rsecOut,whichClusters="mergeClusters"),
               clusterMatrix(mergeOut,whichClusters="mergeClusters"))
  expect_equal(clusterTypes(rsecOut),clusterTypes(mergeOut))
})

test_that("`RSEC` works with hdf5",{
  
  expect_message(rsecOut2<-RSEC(hdfObj, isCount=FALSE,k0s=4:5,
                                reduceMethod="PCA",
                                clusterFunction="hierarchical01", alphas=0.1, nReducedDims=3,
                                subsampleArgs=list(resamp.num=5),
                                seqArgs=list(top.can=0),
                                random.seed=seedValue),
                 "Merging will be done on"
  )
  
  expect_message(rsecOut3<-RSEC(assay(hdfObj), isCount=FALSE,
                                k0s=4:5, reduceMethod="PCA",
                                clusterFunction="hierarchical01", alphas=0.1, nReducedDims=3,
                                subsampleArgs=list(resamp.num=5),
                                seqArgs=list(top.can=0), random.seed=seedValue),
                 "Merging will be done on"
  )
  
  expect_equal(clusterMatrix(rsecOut2),clusterMatrix(rsecOut3))
  
  #no reduce method, do everything on raw data
  #requires numeric/complex matrix/vector arguments
  expect_message(rsecOut1<-RSEC(hdfObj, isCount=FALSE,
                                k0s=4:5,reduceMethod="none",
                                clusterFunction="hierarchical01", alphas=0.1, 
                                seqArgs=list(top.can=0), 
                                subsampleArgs=list(resamp.num=5,clusterFunction="pam"),
                                random.seed=seedValue),
                 "Merging will be done on"
  )
})

