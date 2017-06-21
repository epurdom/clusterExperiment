context("clusterMany")
source("create_objects.R")

test_that("`clusterMany` works with matrix, list of data, ClusterExperiment objects, and
          SummarizedExperiments", {
			  #check all builtin methods
            expect_silent(clustNothing <- clusterMany(mat, ks=c(3,4),clusterFunction=builtInClusterFunctions,
                                       subsample=FALSE, sequential=FALSE,
                                       isCount=FALSE,verbose=FALSE))
			expect_silent(clustDF <- clusterMany(data.frame(mat), ks=c(3,4),clusterFunction=builtInClusterFunctions,
						                                          subsample=FALSE, sequential=FALSE,
						                                          isCount=FALSE,verbose=FALSE))
				   
            expect_is(clustNothing, "ClusterExperiment")
            expect_is(clustNothing, "SummarizedExperiment")
 
            expect_silent(clustNothing2 <- clusterMany(se, ks=c(3,4),clusterFunction="pam",
                                            subsample=FALSE, sequential=FALSE,
                                            isCount=FALSE,verbose=FALSE))
            expect_equal(colData(clustNothing2),colData(se))
            expect_equal(rownames(clustNothing2),rownames(se))
            expect_equal(colnames(clustNothing2),colnames(se))
            expect_equal(metadata(clustNothing2),metadata(se))
            expect_equal(rowData(clustNothing2),rowData(se))

            expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing))
            expect_true(all(clusterTypes(clustNothing)=="clusterMany"))

            #test running on clusterExperiment Object -- should add the new clustering
            expect_silent(clustNothing3 <- clusterMany(ccSE, ks=c(3,4),clusterFunction="pam",
                                         subsample=FALSE, sequential=FALSE,
                                         isCount=FALSE,verbose=FALSE))
            expect_true(nClusters(clustNothing3) == nClusters(ccSE) + 2)
            expect_equal(colData(clustNothing3),colData(ccSE))
            expect_equal(rownames(clustNothing3),rownames(ccSE))
            expect_equal(colnames(clustNothing3),colnames(ccSE))
            expect_equal(metadata(clustNothing3),metadata(ccSE))
            expect_equal(rowData(clustNothing3),rowData(ccSE))
            
            expect_silent(test <- clusterSingle(se,  subsample=FALSE, sequential=FALSE, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=4)),isCount=FALSE))
            expect_silent(clustNothing3<- clusterMany(test, ks=c(3,4),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE))
            expect_silent(clustNothing4<- clusterMany(clustNothing3, ks=c(3:4),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE,eraseOld=TRUE))
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
  expect_silent(cc <- clusterMany(mat, ks=c(3,4),nVarDim=c(10,15),nPCADim=c(3,4),dimReduce=c("none","PCA","var","cv","mad"),clusterFunction="pam",
                    subsample=FALSE, sequential=FALSE,verbose=FALSE,
                    isCount=FALSE)
					)
  #check giving paramMatrix
  expect_silent(param <- clusterMany(mat, ks=c(3,4),nVarDim=c(10,15),nPCADim=c(3,4),dimReduce=c("none","PCA","var"),clusterFunction="pam",
                       subsample=FALSE, sequential=FALSE,run=FALSE,verbose=FALSE,
                       isCount=FALSE))
  #             cc2 <- clusterMany(mat, ks=c(3,4),nVarDim=c(10, 15),nPCADim=c(3,4),dimReduce=c("none","PCA","var"),clusterFunction="pam",
  #                                            subsample=FALSE, sequential=FALSE,verbose=FALSE,
  #                                            isCount=FALSE,paramMatrix=param$paramMatrix,mainClusterArgs=param$mainClusterArgs,seqArgs=param$seqArgs,subsampleArgs=param$subsampleArgs)
  #             expect_equal(cc,cc2)
  
#   #check giving distance -- this still doesn't work. 
#   dist1<-function(x){dist(x,method="manhattan")}
#   dist2<-function(x){dist(x)} ## problem to just give dist because need to grab from global environment
#   cc <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
#                     distFunction=c("dist1","dist2",NA),
#                     subsample=FALSE, sequential=FALSE,verbose=FALSE,
#                     isCount=FALSE)
  
  #check doesn't spit out warnings because alphas/mainClustering args not match 
  expect_silent(clusterMany(mat, clusterFunction=c("pam","hierarchical01"),ks=c(3,4),
                    alphas=c(0.1,0.2),
                    subsample=FALSE, sequential=FALSE,verbose=FALSE,
                    mainClusterArgs=list(clusterArgs=list(evalClusterMethod="average")),
                    isCount=FALSE))
  
  #check doesn't spit out warnings because alphas/mainClustering args not match 
  expect_silent(clusterMany(mat, clusterFunction=c("pam","hierarchical01"),ks=c(3,4),
                            betas=c(.7,.9), minSizes=c(3,5),
                            subsample=FALSE, sequential=FALSE,verbose=FALSE,
                            mainClusterArgs=list(clusterArgs=list(evalClusterMethod="average")),
                            isCount=FALSE))
})

