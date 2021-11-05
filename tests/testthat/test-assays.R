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
    expect_silent(cl1 <- clusterSingle(multi_se, 
      whichAssay = "counts",
      mainClusterArgs=list(clusterArgs=list(k=3),
          clusterFunction="pam"), 
      subsample=FALSE, sequential=FALSE))

    ###Apply to CE
    expect_silent(cl2 <- clusterSingle(multi_cc, 
        whichAssay = "counts",
        mainClusterArgs=list(clusterArgs=list(k=3), 
            clusterFunction="pam"), 
        subsample=FALSE, 
        sequential=FALSE)
    )

    expect_equal(primaryCluster(cl1), primaryCluster(cl2))


    expect_silent(cl3 <- clusterSingle(multi_cc, 
        whichAssay = "logcounts", 
        mainClusterArgs=list(clusterArgs=list(k=3),
            clusterFunction="pam"),
        subsample=FALSE, sequential=FALSE))

    expect_false(all(primaryCluster(cl2) == primaryCluster(cl3)))

    expect_error(clusterSingle(multi_cc, whichAssay = "normalized",
        mainClusterArgs=list(clusterArgs=list(k=3),
            clusterFunction="pam"),
        subsample=FALSE, sequential=FALSE),
        "'normalized' not in names")

    expect_error(clusterSingle(multi_cc, whichAssay = 3,
        mainClusterArgs=list(clusterArgs=list(k=3),
            clusterFunction="pam"),
        subsample=FALSE, sequential=FALSE),
        "subscript is out of bounds")

})

test_that("clusterMany works with non default assays", {
  expect_silent(cl1 <- clusterMany(multi_cc, ks=c(3,4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,
                                   whichAssay = "counts",verbose=FALSE))

  expect_silent(cl2 <- clusterMany(multi_se, ks=c(3,4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,
                                   whichAssay = "counts",verbose=FALSE))

  expect_equal(primaryCluster(cl1), primaryCluster(cl2))


  expect_silent(cl3 <- clusterMany(multi_cc, ks=c(3,4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,
                                   whichAssay = "logcounts",verbose=FALSE))

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

test_that("plotting works wih non default assays", {
  expect_silent(plotFeatureBoxplot(object=multi_cc,feature=1,whichAssay=1))
  expect_silent(plotFeatureBoxplot(object=multi_cc,feature=1,whichAssay=2))

  small_cc <- multi_cc[1:10,c(1:3, 290:293)]
  small_cc <- makeDendrogram(small_cc)
  expect_silent(plotHeatmap(small_cc, whichAssay=1))
  expect_silent(plotHeatmap(small_cc, whichAssay=2))

})


