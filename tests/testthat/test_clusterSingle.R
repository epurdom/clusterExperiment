library(clusterExperiment)

test_that("`clusterSingle` works with matrix, ClusterExperiment objects, and
          SummarizedExperiments", {
            clustNothing <- clusterSingle(simData, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isCount=FALSE)
            expect_is(clustNothing, "ClusterExperiment")
            expect_is(clustNothing, "SummarizedExperiment")


            clustNothing2 <- clusterSingle(se, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isCount=FALSE)
            expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing))

            #test running on clusterExperiment Object -- should add the new clustering
            clustNothing3 <- clusterSingle(clustNothing2, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=4),is=FALSE)
            expect_equal(NCOL(clusterMatrix(clustNothing3)),2)
            expect_equal(length(table(primaryCluster(clustNothing3))),4,info="Check reset primary cluster after run clusterSingle")

          })

test_that("Different options of `clusterSingle` ", {
  #check isCount
  clustCount <- clusterSingle(simCount, clusterFunction="pam",
                           subsample=FALSE, sequential=FALSE,
                           clusterDArgs=list(k=3),isCount=TRUE)
  expect_error(clusterSingle(simData, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,
                          clusterDArgs=list(k=3),isCount=TRUE),info="test error handling for isCount=TRUE when can't take log")

  #check subsample
  clustSubsample <- clusterSingle(simData, clusterFunction="pam",
                               subsample=TRUE, sequential=FALSE,
                               subsampleArgs=list(resamp.num=3),clusterDArgs=list(k=3),isCount=FALSE)
  expect_equal(NCOL(coClustering(clustSubsample)),NCOL(simData))

  #check pca reduction
  clustndims <- clusterSingle(simData, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,dimReduce="PCA",ndims=3,
                          clusterDArgs=list(k=3),isCount=FALSE)
  expect_error(  clusterSingle(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,dimReduce="PCA",ndims=NROW(simData)+1,
                            clusterDArgs=list(k=3),isCount=FALSE))

  #check var reduction
  clustndims <- clusterSingle(simData, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,dimReduce="mostVar",ndims=3,
                          clusterDArgs=list(k=3),isCount=FALSE)
  expect_error(  clusterSingle(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,dimReduce="mostVar",ndims=NROW(simData)+1,
                            clusterDArgs=list(k=3),isCount=FALSE))
  expect_warning(  clusterSingle(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,dimReduce="none",ndims =3,
                            clusterDArgs=list(k=3),isCount=FALSE))

  #check sequential
  clustSeq <- clusterSingle(simData, clusterFunction="pam",
                         subsample=FALSE, sequential=TRUE,
                         isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE))
  expect_error(  clusterSingle(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,
                            isCount=FALSE)) #must specify k0

  #check warning combinations
  expect_error(  clusterSingle(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,seqArgs=list(verbose=FALSE),
                            isCount=FALSE,clusterDArgs=list("typeAlg"=="K")))
  expect_error(  clusterSingle(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,seqArgs=list(verbose=FALSE),
                            isCount=FALSE,clusterDArgs=list("findBestK"==TRUE)))
  expect_error(clusterSingle(simData, clusterFunction="tight",
                          subsample=FALSE, sequential=FALSE,
                          subsampleArgs=list(resamp.num=3),clusterDArgs=list(k=3),isCount=FALSE))
  expect_warning(clusterSingle(simData, clusterFunction="pam",
                            subsample=TRUE, sequential=FALSE,
                            subsampleArgs=list(resamp.num=3),clusterDArgs=list(k=3),isCount=FALSE))
  expect_error(clusterSingle(simData, clusterFunction="pam",
                          subsample=TRUE, sequential=FALSE,
                          subsampleArgs=list(resamp.num=3),isCount=FALSE)) #must give k
  expect_warning(clusterSingle(simData, clusterFunction="tight",
                            subsample=TRUE, sequential=FALSE,
                            subsampleArgs=list(resamp.num=3,k=3),clusterDArgs=list(findBestK=TRUE),isCount=FALSE)) #K argument not match tight
})
