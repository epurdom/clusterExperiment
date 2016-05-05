context("clusterSingle")
source("create_objects.R")

test_that("`clusterSingle` works with matrix, ClusterExperiment objects, and
          SummarizedExperiments", {
            clustNothing <- clusterSingle(mat, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isCount=FALSE)
            expect_is(clustNothing, "ClusterExperiment")
            expect_is(clustNothing, "SummarizedExperiment")


            clustNothing2 <- clusterSingle(se, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isCount=FALSE)
            expect_equal(clusterMatrix(clustNothing2), clusterMatrix(clustNothing))

            #test running on clusterExperiment Object -- should add the new clustering
            clustNothing3 <- clusterSingle(clustNothing2, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=4),is=FALSE)
            expect_equal(NCOL(clusterMatrix(clustNothing3)),2)
            expect_equal(length(table(primaryCluster(clustNothing3))),4,info="Check reset primary cluster after run clusterSingle")

          })

test_that("Different options of `clusterSingle` ", {
  #check isCount
  clustCount <- clusterSingle(smSimCount, clusterFunction="pam",
                           subsample=FALSE, sequential=FALSE,
                           clusterDArgs=list(k=3),isCount=TRUE)
  expect_error(clusterSingle(smSimData, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,
                          clusterDArgs=list(k=3),isCount=TRUE),info="test error handling for isCount=TRUE when can't take log")

  #check subsample
  clustSubsample <- clusterSingle(mat, clusterFunction="pam",
                               subsample=TRUE, sequential=FALSE,
                               subsampleArgs=list(resamp.num=3, k=3),
                               clusterDArgs=list(k=3),isCount=FALSE)
  expect_equal(NCOL(coClustering(clustSubsample)),NCOL(mat))

  #check pca reduction
  clustndims <- clusterSingle(mat, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE, dimReduce="PCA",
                          ndims=3, clusterDArgs=list(k=3),isCount=FALSE)
  expect_error(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE, dimReduce="PCA",
                            ndims=NROW(simData)+1,
                            clusterDArgs=list(k=3),isCount=FALSE))

  #check var reduction
  clustndims <- clusterSingle(mat, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,
                          dimReduce="mostVar", ndims=3,
                          clusterDArgs=list(k=3), isCount=FALSE)
  expect_error(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,
                            dimReduce="mostVar", ndims=NROW(mat)+1,
                            clusterDArgs=list(k=3),isCount=FALSE),
               "the number of most variable features must be strictly less than the number of rows of input data matrix")
  expect_warning(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,
                            dimReduce="none",ndims =3,
                            clusterDArgs=list(k=3),isCount=FALSE),
                 "specifying ndims has no effect if dimReduce==`none`")

  #check sequential
  clustSeq <- clusterSingle(mat, clusterFunction="pam",
                         subsample=FALSE, sequential=TRUE,
                         isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE))
  expect_error(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,
                            isCount=FALSE), "must give seqArgs so as to identify k0")

  #check errors and warnings
  expect_error(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,
                            seqArgs=list(verbose=FALSE),
                            isCount=FALSE,clusterDArgs=list("typeAlg"=="K")),
               "seqArgs must contain element 'k0'")
  expect_error(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=TRUE,
                            seqArgs=list(verbose=FALSE),
                            isCount=FALSE, clusterDArgs=list("findBestK"==TRUE)),
               "seqArgs must contain element 'k0'")
  expect_error(clusterSingle(mat, clusterFunction="tight",
                          subsample=FALSE, sequential=FALSE,
                          subsampleArgs=list(resamp.num=3),
                          clusterDArgs=list(k=3), isCount=FALSE),
               "If not subsampling, clusterFunction must be 'pam'")
  expect_warning(clusterSingle(mat, clusterFunction="pam",
                            subsample=TRUE, sequential=FALSE,
                            subsampleArgs=list(resamp.num=3),
                            clusterDArgs=list(k=3), isCount=FALSE),
                 "did not give 'k' in 'subsampleArgs'.")
  expect_error(clusterSingle(mat, clusterFunction="pam",
                          subsample=TRUE, sequential=FALSE,
                          subsampleArgs=list(resamp.num=3), isCount=FALSE),
               "must pass 'k' in subsampleArgs")
  expect_warning(clusterSingle(mat, clusterFunction="tight",
                            subsample=TRUE, sequential=FALSE,
                            subsampleArgs=list(resamp.num=3,k=3),
                            clusterDArgs=list(findBestK=TRUE),isCount=FALSE),
                 "do not match the choice of typeAlg")
})

test_that("`clusterSingle` preserves the colData and rowData of SE", {
  cl <- clusterSingle(se, clusterFunction="pam",
                      subsample=FALSE, sequential=FALSE,
                      clusterDArgs=list(k=3),isCount=FALSE)

  expect_equal(colData(cl),colData(se))
  expect_equal(rownames(cl),rownames(se))
  expect_equal(colnames(cl),colnames(se))
  expect_equal(metadata(cl),metadata(se))
  expect_equal(rowData(cl),rowData(se))

})
