library(clusterCells)
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version. Can delete this once package is stabilized.
test_that("`compareChoices` works with matrix, list of data, ClusterCells objects, and
          SummarizedExperiments", {
            clustNothing <- compareChoices(simData, ks=c(3,4),clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       isCount=FALSE)
            expect_is(clustNothing, "ClusterCells")
            expect_is(clustNothing, "SummarizedExperiment")
            
            
            se <- SummarizedExperiment(simData)
            clustNothing2 <- compareChoices(se, ks=c(3,4),clusterFunction="pam",
                                            subsample=FALSE, sequential=FALSE,
                                            isCount=FALSE)
            expect_equal(clustNothing, clustNothing2)
            
            #test running on clusterCells Object -- should add the new clustering
            #not yet implemented
            test <- clusterAll(se, clusterFunction="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        clusterDArgs=list(k=4),is=FALSE)
            clustNothing3<- compareChoices(test, ks=c(3,4),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,
                                           isCount=FALSE)
            clustNothing4<- compareChoices(clustNothing3, ks=c(3:4),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,
                                           isCount=FALSE,eraseOld=TRUE) #hitting error. Need to go back (infinite recursion -- why?)
            expect_equal(clustNothing3,clustNothing4)
            clustNothing5<- compareChoices(clustNothing3, ks=c(5:6),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,
                                           isCount=FALSE,eraseOld=FALSE)
            expect_equal(NCOL(allClusters(clustNothing5)),5)
            ppIndex<-pipelineClusterIndex(clustNothing5,print=FALSE)
            expect_equal(as.numeric(table(ppIndex[,"iteration"])),c(2,2))
            #check dim reduce 
            cc <- compareChoices(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,
                                           isCount=FALSE)
            #check giving paramMatrix
            param <- compareChoices(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,run=FALSE,
                                           isCount=FALSE)
            cc2 <- compareChoices(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,
                                           isCount=FALSE,paramMatrix=param)
            expect_equal(cc,cc2)
            
          })
