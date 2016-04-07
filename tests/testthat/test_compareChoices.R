library(clusterCells)
data(simData)

test_that("`compareChoices` works with matrix, list of data, ClusterCells objects, and
          SummarizedExperiments", {
            clustNothing <- compareChoices(simData, ks=c(3,4),clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       isCount=FALSE)
            expect_is(clustNothing, "ClusterCells")
            expect_is(clustNothing, "SummarizedExperiment")
            
            
            se <- SummarizedExperiment(simData)
            clustNothing2 <- clusterAll(se, clusterFunction="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        clusterDArgs=list(k=3),isCount=FALSE)
            expect_equal(clustNothing, clustNothing2)
            
            #test running on clusterCells Object -- should add the new clustering
            clustNothing3 <- clusterAll(clustNothing2, clusterFunction="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        clusterDArgs=list(k=4),is=FALSE)
            expect_equal(NCOL(allClusters(clustNothing3)),2)
            expect_equal(length(table(primaryCluster(clustNothing3))),4,info="Check reset primary cluster after run clusterAll")
            
            })
