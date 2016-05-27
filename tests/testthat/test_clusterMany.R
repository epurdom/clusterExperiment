context("clusterMany")
source("create_objects.R")

test_that("`clusterMany` works with matrix, list of data, ClusterExperiment objects, and
          SummarizedExperiments", {
            clustNothing <- clusterMany(mat, ks=c(3,4),clusterFunction=c("pam","hierarchicalK","hierarchical01","tight"),
                                       subsample=FALSE, sequential=FALSE,
                                       isCount=FALSE,verbose=FALSE)
            expect_is(clustNothing, "ClusterExperiment")
            expect_is(clustNothing, "SummarizedExperiment")
 
            clustNothing2 <- clusterMany(se, ks=c(3,4),clusterFunction="pam",
                                            subsample=FALSE, sequential=FALSE,
                                            isCount=FALSE,verbose=FALSE)
            expect_equal(colData(clustNothing2),colData(se))
            expect_equal(rownames(clustNothing2),rownames(se))
            expect_equal(colnames(clustNothing2),colnames(se))
            expect_equal(metadata(clustNothing2),metadata(se))
            expect_equal(rowData(clustNothing2),rowData(se))

            expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing))
            expect_true(all(clusterTypes(clustNothing)=="clusterMany"))

            #test running on clusterExperiment Object -- should add the new clustering
            clustNothing3 <- clusterMany(ccSE, ks=c(3,4),clusterFunction="pam",
                                         subsample=FALSE, sequential=FALSE,
                                         isCount=FALSE,verbose=FALSE)
            expect_true(nClusters(clustNothing3) == nClusters(ccSE) + 2)
            expect_equal(colData(clustNothing3),colData(ccSE))
            expect_equal(rownames(clustNothing3),rownames(ccSE))
            expect_equal(colnames(clustNothing3),colnames(ccSE))
            expect_equal(metadata(clustNothing3),metadata(ccSE))
            expect_equal(rowData(clustNothing3),rowData(ccSE))
            
            test <- clusterSingle(se, clusterFunction="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        clusterDArgs=list(k=4),isCount=FALSE)
            clustNothing3<- clusterMany(test, ks=c(3,4),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE)
            clustNothing4<- clusterMany(clustNothing3, ks=c(3:4),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE,eraseOld=TRUE)
            expect_equal(clustNothing3,clustNothing4)

            clustNothing5<- clusterMany(clustNothing3, ks=c(5:6),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE,eraseOld=FALSE)
            expect_equal(NCOL(clusterMatrix(clustNothing5)),5)

            ppIndex<-workflowClusterDetails(clustNothing5)
            expect_equal(as.numeric(table(ppIndex[,"iteration"])),c(2,2))

 
          })
test_that("`clusterMany` works changing parameters", {
  #check dim reduce
  cc <- clusterMany(mat, ks=c(3,4),nVarDim=c(10,15),nPCADim=c(3,4),dimReduce=c("none","PCA","var","cv","mad"),clusterFunction="pam",
                    subsample=FALSE, sequential=FALSE,verbose=FALSE,
                    isCount=FALSE)
  #check giving paramMatrix
  param <- clusterMany(mat, ks=c(3,4),nVarDim=c(10,15),nPCADim=c(3,4),dimReduce=c("none","PCA","var"),clusterFunction="pam",
                       subsample=FALSE, sequential=FALSE,run=FALSE,verbose=FALSE,
                       isCount=FALSE)
  #             cc2 <- clusterMany(mat, ks=c(3,4),nVarDim=c(10, 15),nPCADim=c(3,4),dimReduce=c("none","PCA","var"),clusterFunction="pam",
  #                                            subsample=FALSE, sequential=FALSE,verbose=FALSE,
  #                                            isCount=FALSE,paramMatrix=param$paramMatrix,clusterDArgs=param$clusterDArgs,seqArgs=param$seqArgs,subsampleArgs=param$subsampleArgs)
  #             expect_equal(cc,cc2)
  
#   #check giving distance -- this still doesn't work. 
#   dist1<-function(x){dist(x,method="manhattan")}
#   dist2<-function(x){dist(x)} ## problem to just give dist because need to grab from global environment
#   cc <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
#                     distFunction=c("dist1","dist2",NA),
#                     subsample=FALSE, sequential=FALSE,verbose=FALSE,
#                     isCount=FALSE)
  
  #check doesn't spit out warnings because alphas/clusterD args not match 
  expect_silent(clusterMany(mat, clusterFunction=c("pam","hierarchical01"),ks=c(3,4),
                    alphas=c(0.1,0.2),
                    subsample=FALSE, sequential=FALSE,verbose=FALSE,
                    clusterDArgs=list(clusterArgs=list(evalClusterMethod="average")),
                    isCount=FALSE))
})
