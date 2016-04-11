library(clusterCells)
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version. Can delete this once package is stabilized.
test_that("`clusterMany` works with matrix, list of data, ClusterCells objects, and
          SummarizedExperiments", {
            clustNothing <- clusterMany(simData, ks=c(3,4),clusterMethod="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       isCount=FALSE,verbose=FALSE)
            expect_is(clustNothing, "ClusterCells")
            expect_is(clustNothing, "SummarizedExperiment")
            clusterLabels(clustNothing,whichClusters="pipeline")
            
            se <- SummarizedExperiment(simData)
            clustNothing2 <- clusterMany(se, ks=c(3,4),clusterMethod="pam",
                                            subsample=FALSE, sequential=FALSE,
                                            isCount=FALSE,verbose=FALSE)
            expect_equal(clustNothing, clustNothing2)
            
            #test running on clusterCells Object -- should add the new clustering
            #not yet implemented
            test <- clusterAll(se, clusterFunction="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        clusterDArgs=list(k=4),isCount=FALSE)
            clustNothing3<- clusterMany(test, ks=c(3,4),clusterMethod="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE)
            clustNothing4<- clusterMany(clustNothing3, ks=c(3:4),clusterMethod="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE,eraseOld=TRUE) 
            expect_equal(clustNothing3,clustNothing4)
            clustNothing5<- clusterMany(clustNothing3, ks=c(5:6),clusterMethod="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE,eraseOld=FALSE)
            expect_equal(NCOL(allClusters(clustNothing5)),5)
            ppIndex<-pipelineClusterDetails(clustNothing5)
            expect_equal(as.numeric(table(ppIndex[,"iteration"])),c(2,2))
            #check dim reduce 
            cc <- clusterMany(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterMethod="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE)
            #check giving paramMatrix
            param <- clusterMany(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterMethod="pam",
                                           subsample=FALSE, sequential=FALSE,run=FALSE,verbose=FALSE,
                                           isCount=FALSE)
            cc2 <- clusterMany(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterMethod="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE,paramMatrix=param$paramMatrix,clusterDArgs=param$clusterDArgs,seqArgs=param$seqArgs,subsampleArgs=param$subsampleArgs)
            expect_equal(cc,cc2)
            
          })
