library(clusterExperiment)
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version. Can delete this once package is stabilized.

test_that("`combineMany` works with matrix and ClusterExperiment objects", {
            clustNothing <- clusterMany(simData, ks=c(3,4),clusterFunction="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        isCount=FALSE,verbose=FALSE)

            expect_error(combineMany(clustNothing), "missing")

            shared1 <- combineMany(clusterMatrix(clustNothing))
            shared2 <- combineMany(clustNothing, "all")
            expect_equal(shared1$clustering, primaryCluster(shared2))

            shared3 <- combineMany(clustNothing, "pipeline")
            expect_equal(shared2, shared3)

            shared4 <- combineMany(clustNothing, 1:nClusters(clustNothing))
            expect_equal(shared3, shared4)

            shared5 <- combineMany(clustNothing, "pipeline",
                                          proportion=.5)

            shared6 <- combineMany(clustNothing, "pipeline",
                                          proportion=.5,
                                          clusterFunction="tight")

            expect_error(combineMany(clustNothing, "pipeline",
                                          proportion=.5,
                                          clusterFunction="pam"), "implemented")

            expect_true("combineMany" %in% clusterType(shared6))
            expect_true("combineMany" %in% clusterLabels(shared6))
            expect_true(all(primaryCluster(shared6)==clusterMatrix(shared6)[,"combineMany"]))
})

test_that("`combineMany` works when multiple runs of pipeline", {
  clustNothing <- clusterMany(simData, ks=c(3,4),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE)

  shared1 <- combineMany(clustNothing, "all")

  clustNothing2 <- clusterMany(shared1, ks=c(5,6), clusterFunction="pam",
                               subsample=FALSE, sequential=FALSE,
                               isCount=FALSE,verbose=FALSE)

  shared3 <- combineMany(clustNothing2, "all")
  shared4 <- combineMany(clusterMatrix(clustNothing2))
  expect_equal(shared4$clustering, primaryCluster(shared3))

  shared5 <- combineMany(clustNothing2, "pipeline")
  shared6 <- combineMany(clusterMatrix(clustNothing2)[,1:2])
  expect_equal(shared6$clustering, primaryCluster(shared5))


  clustNothing3 <- addClusters(clustNothing2, primaryCluster(shared5))
  shared7 <- combineMany(clustNothing3, "all")
  shared8 <- combineMany(clustNothing3, "pipeline")
})
