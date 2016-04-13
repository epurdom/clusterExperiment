library(clusterExperiment)
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version. Can delete this once package is stabilized.

test_that("`findSharedClusters` works with matrix and ClusterExperiment objects", {
            clustNothing <- clusterMany(simData, ks=c(3,4),clusterMethod="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        isCount=FALSE,verbose=FALSE)

            expect_error(findSharedClusters(clustNothing), "missing")

            shared1 <- findSharedClusters(allClusters(clustNothing))
            shared2 <- findSharedClusters(clustNothing, "all")
            expect_equal(shared1$clustering, primaryCluster(shared2))

            shared3 <- findSharedClusters(clustNothing, "pipeline")
            expect_equal(shared2, shared3)

            shared4 <- findSharedClusters(clustNothing, 1:nClusters(clustNothing))
            expect_equal(shared3, shared4)

            shared5 <- findSharedClusters(clustNothing, "pipeline",
                                          proportion=.5)

            shared6 <- findSharedClusters(clustNothing, "pipeline",
                                          proportion=.5,
                                          clusterFunction="tight")

            expect_error(findSharedClusters(clustNothing, "pipeline",
                                          proportion=.5,
                                          clusterFunction="pam"), "implemented")

            expect_true("findSharedClusters" %in% clusterType(shared6))
            expect_true("findSharedClusters" %in% clusterLabels(shared6))
            expect_true(all(primaryCluster(shared6)==allClusters(shared6)[,"findSharedClusters"]))
})

test_that("`findSharedClusters` works when multiple runs of pipeline", {
  clustNothing <- clusterMany(simData, ks=c(3,4),clusterMethod="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE)

  shared1 <- findSharedClusters(clustNothing, "all")

  clustNothing2 <- clusterMany(shared1, ks=c(5,6), clusterMethod="pam",
                               subsample=FALSE, sequential=FALSE,
                               isCount=FALSE,verbose=FALSE)

  shared3 <- findSharedClusters(clustNothing2, "all")
  shared4 <- findSharedClusters(allClusters(clustNothing2))
  expect_equal(shared4$clustering, primaryCluster(shared3))

  shared5 <- findSharedClusters(clustNothing2, "pipeline")
  shared6 <- findSharedClusters(allClusters(clustNothing2)[,1:2])
  expect_equal(shared6$clustering, primaryCluster(shared5))


  clustNothing3 <- addClusters(clustNothing2, primaryCluster(shared5))
  shared7 <- findSharedClusters(clustNothing3, "all")
  shared8 <- findSharedClusters(clustNothing3, "pipeline")
})
