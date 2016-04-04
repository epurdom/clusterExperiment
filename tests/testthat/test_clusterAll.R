library(clusterCells)
data(simData)

test_that("`clusterAll` works with matrix, ClusterCells objects, and
          SummarizedExperiments", {
            clustNothing <- clusterAll(simData, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isCount=FALSE)
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
            

  #check isCount
            clustCount <- clusterAll(simCount, clusterFunction="pam",
                                     subsample=FALSE, sequential=FALSE,
                                     clusterDArgs=list(k=3),isCount=TRUE)
            expect_error(clusterAll(simData, clusterFunction="pam",
                                    subsample=FALSE, sequential=FALSE,
                                    clusterDArgs=list(k=3),isCount=TRUE),info="test error handling for isCount=TRUE when can't take log")

            #check subsample
            clustSubsample <- clusterAll(simData, clusterFunction="pam",
                                     subsample=TRUE, sequential=FALSE,
                                     subsampleArgs=list(resamp.num=3),clusterDArgs=list(k=3),isCount=FALSE)
            expect_equal(NCOL(coClustering(clustSubsample)),NCOL(simData))
            #check npcs
            clustNPCS <- clusterAll(simData, clusterFunction="pam",
                                         subsample=FALSE, sequential=FALSE,npcs=3,
                                         clusterDArgs=list(k=3),isCount=FALSE)
           expect_error(  clusterAll(simData, clusterFunction="pam",
                                                  subsample=FALSE, sequential=FALSE,npcs=NROW(simData)+1,
                                                  clusterDArgs=list(k=3),isCount=FALSE))
            
            
            
          })
