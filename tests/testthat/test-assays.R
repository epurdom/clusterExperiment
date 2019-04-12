context("Assays")

## Construct ClusterExperiment with multiple named assays
multi_se <- SummarizedExperiment(assays = list(counts = simCount,
                                               logcounts = log1p(simCount)))
multi_cc <- ClusterExperiment(multi_se, trueCluster)
ccTrue <- ClusterExperiment(simCount, trueCluster)
ccTrue2 <- ClusterExperiment(simCount, trueCluster, transformation = log1p)
seedValue<-495 #01875 works for sample.kind="Rejection"

test_that("clusterSingle works with non defalt assays", {

  ###Apply to SE
  expect_silent(cl1 <- clusterSingle(multi_se, whichAssay = "counts",
                                     mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"), subsample=FALSE, sequential=FALSE))

  ###Apply to CE
  expect_silent(cl2 <- clusterSingle(multi_cc, whichAssay = "counts",
                                     mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"), subsample=FALSE, sequential=FALSE))

  expect_equal(primaryCluster(cl1), primaryCluster(cl2))


  expect_silent(cl3 <- clusterSingle(multi_cc, whichAssay = "logcounts",
                                     mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
                                     subsample=FALSE, sequential=FALSE))

  expect_false(all(primaryCluster(cl2) == primaryCluster(cl3)))

  expect_error(clusterSingle(multi_cc, whichAssay = "normalized",
                             mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
                             subsample=FALSE, sequential=FALSE),
               "'normalized' not in names")

  expect_error(clusterSingle(multi_cc, whichAssay = 3,
                             mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
                             subsample=FALSE, sequential=FALSE),
               "subscript is out of bounds")

})

test_that("clusterMany works with non default assays", {

  expect_silent(cl1 <- clusterMany(multi_cc, ks=c(3,4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,
                                   whichAssay = "counts"))

  expect_silent(cl2 <- clusterMany(multi_se, ks=c(3,4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,
                                   whichAssay = "counts"))

  expect_equal(primaryCluster(cl1), primaryCluster(cl2))


  expect_silent(cl3 <- clusterMany(multi_cc, ks=c(3,4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,
                                   whichAssay = "logcounts"))

  expect_false(all(primaryCluster(cl1) == primaryCluster(cl3)))
})

test_that("mergeClusters works with non default assays", {

  expect_silent(cl2 <- clusterSingle(multi_cc, whichAssay = "counts",
                                     mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
                                     subsample=FALSE, sequential=FALSE))

  expect_silent(cl3 <- clusterSingle(multi_cc, whichAssay = "logcounts",
                                     mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
                                     subsample=FALSE, sequential=FALSE))

  expect_silent(cl3 <- makeDendrogram(cl3, whichAssay = "logcounts"))
  expect_message(merged <- mergeClusters(x=cl3, whichAssay = "logcounts",DEMethod="limma"),
                 "Merging will be done on")

})

test_that("RSEC works wih non default assays", {
  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "counts"),
                 "Merging will be done on")

  expect_message(out2<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")

  expect_message(out3<-RSEC(x=ccTrue, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue),
                 "Merging will be done on")

  expect_message(out4<-RSEC(x=ccTrue2, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue),
                 "Merging will be done on")

  expect_false(all(primaryCluster(out1) == primaryCluster(out2)))

  expect_equivalent(out1, out3)
  expect_equivalent(out2, out4)
})

test_that("plotting works wih non default assays", {
  expect_silent(plotFeatureBoxplot(object=multi_cc,feature=1,whichAssay=1))
  expect_silent(plotFeatureBoxplot(object=multi_cc,feature=1,whichAssay=2))

  small_cc <- multi_cc[1:10,c(1:3, 290:293)]
  small_cc <- makeDendrogram(small_cc)
  expect_silent(plotHeatmap(small_cc, whichAssay=1))
  expect_silent(plotHeatmap(small_cc, whichAssay=2))

})

test_that("RSEC works independent of assay order", {
  multi_se <- SummarizedExperiment(assays = list(counts = simCount,
                                                 logcounts = log1p(simCount)))
  multi_se2 <- SummarizedExperiment(assays = list(logcounts = log1p(simCount),
                                                  counts = simCount))
  multi_cc <- ClusterExperiment(multi_se, trueCluster)
  multi_cc2 <- ClusterExperiment(multi_se2, trueCluster)

  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")

  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")

  expect_equal(out1, out2)

  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = 2),
                 "Merging will be done on")

  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="none",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = 1),
                 "Merging will be done on")

  expect_equal(out1, out2)

  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="PCA", nReducedDims = 50,
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "counts"),
                 "Merging will be done on")

  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="PCA", nReducedDims = 50,
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "counts"),
                 "Merging will be done on")

  expect_equal(out1, out2)

  expect_message(out1<-RSEC(x=multi_cc, reduceMethod="var",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")

  expect_message(out2<-RSEC(x=multi_cc2, reduceMethod="var",
                            k0s=4:5, clusterFunction="tight", alphas=0.1,
                            betas=0.9, dendroReduce="none", minSizes=1,
                            subsampleArgs=list(resamp.num=5),
                            random.seed=seedValue, whichAssay = "logcounts"),
                 "Merging will be done on")

  expect_equal(out1, out2)

})
