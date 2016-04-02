library(clusterCells)
data(simData)

test_that("`clusterAll` works with matrix, ClusterCells objects, and
          SummarizedExperiments", {
            clustNothing <- clusterAll(simData, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isLog=TRUE)
            se <- SummarizedExperiment(simData)
            clustNothing2 <- clusterAll(se, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isLog=TRUE)
            clustNothing3 <- clusterAll(clustNothing2, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3))

            expect_equal(clustNothing, clustNothing2)
            expect_equal(clustNothing, clustNothing3)

            expect_is(clustNothing, "ClusterCells")
            expect_is(clustNothing, "SummarizedExperiment")
          })
