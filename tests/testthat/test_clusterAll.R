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
            
          })

test_that("Different options of `clusterAll` ", {
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
  
  #check pca reduction
  clustndims <- clusterAll(simData, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,dimReduce="PCA",ndims=3,
                          clusterDArgs=list(k=3),isCount=FALSE)
  expect_error(  clusterAll(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,dimReduce="PCA",ndims=NROW(simData)+1,
                            clusterDArgs=list(k=3),isCount=FALSE))
  
  #check var reduction
  clustndims <- clusterAll(simData, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,dimReduce="mostVar",ndims=3,
                          clusterDArgs=list(k=3),isCount=FALSE)
  expect_error(  clusterAll(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,dimReduce="mostVar",ndims=NROW(simData)+1,
                            clusterDArgs=list(k=3),isCount=FALSE))
  expect_warning(  clusterAll(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,dimReduce="none",ndims =3,
                            clusterDArgs=list(k=3),isCount=FALSE))
  
  #check sequential
  clustSeq <- clusterAll(simData, clusterFunction="pam",
                         subsample=FALSE, sequential=TRUE,ndims=3,
                         isCount=FALSE,seqArgs=list(k0=5))
  expect_error(  clusterAll(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,ndims=3,
                            isCount=FALSE)) #must specify k0

  #check warning combinations
  expect_error(  clusterAll(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,
                            isCount=FALSE,clusterDArgs=list("typeAlg"=="K"))) 
  expect_error(  clusterAll(simData, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,
                            isCount=FALSE,clusterDArgs=list("findBestK"==TRUE))) 
  expect_error(clusterAll(simData, clusterFunction="tight",
                          subsample=FALSE, sequential=FALSE,
                          subsampleArgs=list(resamp.num=3),clusterDArgs=list(k=3),isCount=FALSE))
  expect_warning(clusterAll(simData, clusterFunction="pam",
                            subsample=TRUE, sequential=FALSE,
                            subsampleArgs=list(resamp.num=3),clusterDArgs=list(k=3),isCount=FALSE))
  expect_error(clusterAll(simData, clusterFunction="pam",
                          subsample=TRUE, sequential=FALSE,
                          subsampleArgs=list(resamp.num=3),isCount=FALSE)) #must give k
  expect_warning(clusterAll(simData, clusterFunction="tight",
                            subsample=TRUE, sequential=FALSE,
                            subsampleArgs=list(resamp.num=3,k=3),clusterDArgs=list(findBestK=TRUE),isCount=FALSE)) #K argument not match tight
})