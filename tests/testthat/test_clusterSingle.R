context("clusterSingle")
source("create_objects.R")


test_that("`clusterSingle` works with matrix, ClusterExperiment objects, and
          SummarizedExperiments", {
            clustNothing <- clusterSingle(mat, clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(k=3),isCount=FALSE)
            expect_equal(clusterLabels(clustNothing),"clusterSingle")
            expect_is(clustNothing, "ClusterExperiment")
            expect_is(clustNothing, "SummarizedExperiment")

            #test clusterLabel
            clustNothing2 <- clusterSingle(mat, clusterFunction="pam",
                                          subsample=FALSE, sequential=FALSE,
                                          clusterDArgs=list(k=3),isCount=FALSE,clusterLabel="myownClustering")
            expect_equal(clusterLabels(clustNothing2),"myownClustering")
            
            
            #test default 01 distance
            x1 <- clusterSingle(mat, clusterFunction="tight",
                                          subsample=FALSE, sequential=FALSE,
                                          isCount=FALSE)
            expect_error(clusterSingle(mat, clusterFunction="tight",
                                       subsample=FALSE, sequential=FALSE,
                                       clusterDArgs=list(distFunction=function(x){dist(x,method="manhattan")}),isCount=FALSE),"distance function must give values between 0 and 1")
            
            #test default 01 distance
            x2<-clusterSingle(mat, clusterFunction="tight",
                                       subsample=FALSE, sequential=FALSE,
                                       isCount=FALSE)
            
            #warn wrong arguments
            expect_warning(clusterSingle(mat, clusterFunction="tight",
                              subsample=FALSE, sequential=FALSE,
                              clusterDArgs=list(k=3),isCount=FALSE),"do not match the choice of typeAlg")
            #turn off warning
            expect_silent(clusterSingle(mat, clusterFunction="tight",
                                        subsample=FALSE, sequential=FALSE,
                                        clusterDArgs=list(k=3,checkArgs=FALSE),isCount=FALSE))
            
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

test_that("Different options algorithms of `clusterD` ", {
  #check algorithms
  clusterSingle(mat, clusterFunction="tight",
                subsample=FALSE, sequential=FALSE,
                isCount=FALSE)
  clusterSingle(mat, clusterFunction="hierarchical01",
                subsample=FALSE, sequential=FALSE,
                isCount=FALSE)
  clusterSingle(mat, clusterFunction="hierarchicalK", clusterDArgs=list("k"=3),
                subsample=FALSE, sequential=FALSE,
                isCount=FALSE)
  #K algorithm options
  clusterSingle(mat, clusterFunction="hierarchicalK",
                subsample=FALSE, sequential=FALSE, clusterDArgs=list(findBestK=TRUE,removeSil=TRUE), 
                isCount=FALSE)
  clusterSingle(mat, clusterFunction="pam", clusterDArgs=list(findBestK=TRUE,removeSil=TRUE),
                subsample=FALSE, sequential=FALSE,
                isCount=FALSE)
  
  ########
  #Check clusterD
  ########
  ###Check pam exactly same:
  x<-clusterD(mat, clusterFunction="pam",k=3,
           minSize=1, removeSil=FALSE)
  expect_equal(length(x),ncol(mat))
  x2<-cluster::pam(t(mat),k=3,cluster.only=TRUE)
  expect_equal(x,x2)
  ###Check hierarchicalK exactly same:
  x<-clusterD(mat, clusterFunction="hierarchicalK",k=3,
              minSize=1, removeSil=FALSE)
  expect_equal(length(x),ncol(mat))
  x2<-stats::cutree(stats::hclust(dist(t(mat))),k=3)
  expect_equal(x,x2)
  
  
  #check giving wrong parameters gives warning:
  expect_warning(clusterD(mat, clusterFunction="tight", alpha=0.1,
      minSize=5, removeSil=TRUE),"do not match the choice of typeAlg")
  expect_warning(clusterD(mat, clusterFunction="pam", alpha=0.1,
        minSize=5, removeSil=TRUE, findBestK=TRUE),"do not match the choice of typeAlg")
  expect_warning(clusterD(mat, clusterFunction="tight", alpha=0.1,
    clusterArgs=list(evalClusterMethod="average")),"arguments passed via clusterArgs")
  expect_warning(clusterD(mat, clusterFunction="hierarchical01", alpha=0.1,
   clusterArgs=list(minSize.core=4)),"arguments passed via clusterArgs")
  #check turn off if checkArgs=TRUE
  expect_silent(clusterD(mat, clusterFunction="tight", alpha=0.1,checkArgs=FALSE,
                          minSize=5, removeSil=TRUE))
  expect_silent(clusterD(mat, clusterFunction="pam", alpha=0.1,checkArgs=FALSE,
                          minSize=5, removeSil=TRUE, findBestK=TRUE))
  expect_silent(clusterD(mat, clusterFunction="tight", alpha=0.1,checkArgs=FALSE,
                          clusterArgs=list(evalClusterMethod="average")))
  expect_silent(clusterD(mat, clusterFunction="hierarchical01", alpha=0.1,checkArgs=FALSE,
                          clusterArgs=list(minSize.core=4)))
  
})

test_that("Different options of subsampling",{
    #check subsample
    clustSubsample <- clusterSingle(mat, clusterFunction="pam",
                                    subsample=TRUE, sequential=FALSE,
                                    subsampleArgs=list(resamp.num=3, k=3),
                                    clusterDArgs=list(k=3),isCount=FALSE)
    expect_equal(NCOL(coClustering(clustSubsample)),NCOL(mat))
    clusterSingle(mat, clusterFunction="pam",
                  subsample=TRUE, sequential=FALSE,
                  subsampleArgs=list(resamp.num=3, k=3,clusterFunction="kmeans"),
                  clusterDArgs=list(k=3),isCount=FALSE)
    set.seed(1045)
    clusterSingle(mat, clusterFunction="pam",
                  subsample=TRUE, sequential=FALSE,
                  subsampleArgs=list(resamp.num=20, k=3,classifyMethod="InSample"),
                  clusterDArgs=list(k=3),isCount=FALSE)
    set.seed(1045)
    clusterSingle(mat, clusterFunction="pam",
                  subsample=TRUE, sequential=FALSE,
                  subsampleArgs=list(resamp.num=40, k=3,classifyMethod="OutOfSample"),
                  clusterDArgs=list(k=3),isCount=FALSE)
    set.seed(1045)
    expect_error(clusterSingle(mat, clusterFunction="pam",
                               subsample=TRUE, sequential=FALSE,
                               subsampleArgs=list(resamp.num=20, k=3,classifyMethod="OutOfSample"),
                               clusterDArgs=list(k=3),isCount=FALSE),"NA values found in Dbar")
    
    #errors in missing args in subsample
    expect_warning(clusterSingle(mat, clusterFunction="pam",
                                 subsample=TRUE, sequential=FALSE,
                                 subsampleArgs=list(resamp.num=3),
                                 clusterDArgs=list(k=3), isCount=FALSE),
                   "did not give 'k' in 'subsampleArgs'.")
    expect_error(clusterSingle(mat, clusterFunction="pam",
                               subsample=TRUE, sequential=FALSE,
                               subsampleArgs=list(resamp.num=3), isCount=FALSE),
                 "must pass 'k' in subsampleArgs")
    
})

test_that("Different options of clusterD",{
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
    expect_warning(clusterSingle(mat, clusterFunction="tight",
                                 subsample=FALSE, sequential=FALSE,
                                 clusterDArgs=list(k=3), isCount=FALSE),
                   "do not match the choice of typeAlg")
    expect_warning(clusterSingle(mat, clusterFunction="tight",
                                 subsample=FALSE, sequential=FALSE,
                                 clusterDArgs=list(findBestK=TRUE),isCount=FALSE),
                   "do not match the choice of typeAlg")
})

test_that("Different options of seqCluster",{
    #check sequential
    clustSeq <- clusterSingle(mat, clusterFunction="pam",
                              subsample=FALSE, sequential=TRUE,
                              isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE))
    expect_error(clusterSingle(mat, clusterFunction="pam",
                               subsample=FALSE, sequential=TRUE,
                               isCount=FALSE), "must give seqArgs so as to identify k0")
    
    clustSeq <- clusterSingle(mat, clusterFunction="tight",
                              subsample=FALSE, sequential=TRUE,
                              isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE))
    clustSeq <- clusterSingle(mat, clusterFunction="hierarchicalK",
                              subsample=FALSE, sequential=TRUE,
                              isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE))
    clustSeq <- clusterSingle(mat, clusterFunction="hierarchical01",
                              subsample=FALSE, sequential=TRUE,
                              isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE))
    
    
})

test_that("Different options of `clusterSingle` ", {
  #check isCount
  clustCount <- clusterSingle(smSimCount, clusterFunction="pam",
                           subsample=FALSE, sequential=FALSE,
                           clusterDArgs=list(k=3),isCount=TRUE)
  expect_error(clusterSingle(smSimData, clusterFunction="pam",
                          subsample=FALSE, sequential=FALSE,
                          clusterDArgs=list(k=3),isCount=TRUE),info="test error handling for isCount=TRUE when can't take log")


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
                          dimReduce="var", ndims=3,
                          clusterDArgs=list(k=3), isCount=FALSE)
  expect_error(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,
                            dimReduce="var", ndims=NROW(mat)+1,
                            clusterDArgs=list(k=3),isCount=FALSE),
               "the number of most variable features must be strictly less than the number of rows of input data matrix")
  expect_warning(clusterSingle(mat, clusterFunction="pam",
                            subsample=FALSE, sequential=FALSE,
                            dimReduce="none",ndims =3,
                            clusterDArgs=list(k=3),isCount=FALSE),
                 "specifying ndims has no effect if dimReduce==`none`")

  clustndims <- clusterSingle(mat, clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE, dimReduce="cv",
                              ndims=3, clusterDArgs=list(k=3),isCount=FALSE)
  clustndims <- clusterSingle(mat, clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE, dimReduce="mad",
                              ndims=3, clusterDArgs=list(k=3),isCount=FALSE)
  

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

