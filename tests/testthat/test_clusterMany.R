library(clusterExperiment)
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version. Can delete this once package is stabilized.
sData<-data.frame(sample(letters[2:5],size=NCOL(simData),replace=TRUE),sample(2:5,size=NCOL(simData),replace=TRUE))
sData<-data.frame(sData,sample(LETTERS[2:5],size=NCOL(simData),replace=TRUE),stringsAsFactors=FALSE)
colnames(sData)<-c("A","B","C")

test_that("`clusterMany` works with matrix, list of data, ClusterExperiment objects, and
          SummarizedExperiments", {
            clustNothing <- clusterMany(simData, ks=c(3,4),clusterFunction="pam",
                                       subsample=FALSE, sequential=FALSE,
                                       isCount=FALSE,verbose=FALSE)
            expect_is(clustNothing, "ClusterExperiment")
            expect_is(clustNothing, "SummarizedExperiment")
            clusterLabels(clustNothing,whichClusters="pipeline")

            ##Make this better here:
            se <- SummarizedExperiment(simData,colData=sData)
            clustNothing2 <- clusterMany(se, ks=c(3,4),clusterFunction="pam",
                                            subsample=FALSE, sequential=FALSE,
                                            isCount=FALSE,verbose=FALSE)
            expect_equal(colData(clustNothing2),colData(se)) 
            expect_equal(rownames(clustNothing2),rownames(se)) 
            expect_equal(colnames(clustNothing2),colnames(se)) 
            expect_equal(metadata(clustNothing2),metadata(se)) 
            expect_equal(rowData(clustNothing2),rowData(se)) 
            
            expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing))

            #test running on clusterExperiment Object -- should add the new clustering
            #not yet implemented
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
            ppIndex<-pipelineClusterDetails(clustNothing5)
            expect_equal(as.numeric(table(ppIndex[,"iteration"])),c(2,2))
            #check dim reduce
            cc <- clusterMany(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE)
            #check giving paramMatrix
            param <- clusterMany(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,run=FALSE,verbose=FALSE,
                                           isCount=FALSE)
            cc2 <- clusterMany(simData, ks=c(3,4),nVarDim=c(15,20),nPCADim=c(3,4),dimReduce=c("none","PCA","mostVar"),clusterFunction="pam",
                                           subsample=FALSE, sequential=FALSE,verbose=FALSE,
                                           isCount=FALSE,paramMatrix=param$paramMatrix,clusterDArgs=param$clusterDArgs,seqArgs=param$seqArgs,subsampleArgs=param$subsampleArgs)
            expect_equal(cc,cc2)

          })
