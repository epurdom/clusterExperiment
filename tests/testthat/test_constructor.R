library(clusterExperiment)

mat <- matrix(data=rnorm(200), ncol=10)
mat[1,1]<- -1 #force a negative value
labels <- as.character(gl(5, 2))
labels[c(1:2)]<- c("-1","-2") #make sure some not assigned
labels<-factor(labels)
chLabels<-rep(LETTERS[1:5],each=2)
chLabels[c(2:3)]<- c("-1","-2") #make sure some not assigned
labMat<-cbind(as.numeric(as.character(labels)),as.numeric(as.character(labels)))
se <- SummarizedExperiment(mat)

cc <- clusterExperiment(mat, as.numeric(as.character(labels))+2, transformation = function(x){x})
cc2 <- clusterExperiment(se, as.numeric(as.character(labels)), transformation = function(x){x})
test_that("`clusterExperiment` constructor works with matrix and
          SummarizedExperiments", {
            expect_error(clusterExperiment(mat), "missing")
            expect_error(clusterExperiment(mat,as.numeric(labels),transformation=log),info="Error checking transFun")
            expect_error(clusterExperiment(mat, as.numeric(labels)), "missing")
            expect_error(clusterExperiment(mat, labels[1:2], function(x){x}),
                         "must be a matrix of rows equal")
            expect_error(clusterExperiment(as.data.frame(mat), labels, function(x){x}),
                         "unable to find an inherited method for function")
            #test character input
            ccChar<-clusterExperiment(mat, chLabels, function(x){x})
            expect_is(primaryCluster(ccChar),"numeric")
            expect_is(primaryClusterNamed(ccChar),"character")
            expect_equal(sort(unique(primaryClusterNamed(ccChar))),sort(unique(chLabels)))
            
            #test factor input
            clusterExperiment(mat, labels, function(x){x})

            expect_is(cc, "ClusterExperiment")
            expect_is(cc, "SummarizedExperiment")
            
            expect_equal(nSamples(cc),ncol(mat))
            expect_equal(nFeatures(cc),nrow(mat))
            expect_equal(nClusters(cc),1)

            clusterExperiment(se,labMat,transformation=function(x){x})

                      })

test_that("adding clusters, setting primary labels and remove unclustered cells
          work as promised", {
            expect_equal(NCOL(clusterMatrix(cc)), 1)
            expect_is(transformation(cc),"function")

            c1 <- addClusters(cc, rep(c(-1, 1), each=5))
            expect_equal(NCOL(clusterMatrix(c1)), 2)
            expect_equal(length(clusterType(c1)), 2)
            expect_equal(length(clusterInfo(c1)), 2)
            expect_equal(primaryCluster(c1), primaryCluster(cc))
            primaryClusterIndex(c1) <- 2
            expect_false(all(primaryCluster(c1)==primaryCluster(cc)))

              #check adding a clusterExperimentObject
            c3<-addClusters(cc,cc)
            expect_equal(NCOL(clusterMatrix(c3)), 2)
            expect_equal(length(clusterType(c3)), 2)
            expect_equal(length(clusterInfo(c3)), 2)
            expect_equal(primaryCluster(c3), primaryCluster(cc))
            
            expect_error(clusterLabels(c3)<-c("User","User"))
            clusterLabels(c3)<-c("User1","User2")
            clusterLabels(c3)[1]<-"User4"
            expect_error(clusterLabels(c3)[1]<-"User2","duplicated clusterLabels")
            expect_equal(length(clusterLabels(c3)),2)
            expect_equal(length(clusterLabels(c3,"pipeline")),0) #nothing get back
            
              #check adding matrix of clusters
            c4<-addClusters(cc,cbind(rep(c(-1, 1), each=5),rep(c(2, 1), each=5)),type="New")
            expect_equal(NCOL(clusterMatrix(c4)), 3)
            expect_equal(length(clusterType(c4)), 3)
            expect_equal(length(clusterInfo(c4)), 3)
            expect_equal(primaryCluster(c4), primaryCluster(cc))
            
            ccR<-cc
            dendrogram(ccR)<-NULL
            expect_equal(ccR, removeUnclustered(cc))
#Need to fix:
#             Error: Test failed: 'adding clusters, setting primary labels and remove unclustered cells
#           work as promised'
#             Not expected: ccR not equal to removeUnclustered(cc)
#             Attributes: < Component “clusterLegend”: names for target but not for current >.
#             
            c2 <- removeUnclustered(c1)
            expect_equal(NCOL(c2), 5)
            
            #check removing index of clusters
            c5<-removeClusters(c4,1)
            expect_equal(NCOL(clusterMatrix(c5)), 2)
            expect_equal(length(clusterType(c5)), 2)
            expect_equal(length(clusterInfo(c5)), 2)
            expect_equal(primaryCluster(c4), primaryCluster(removeClusters(c4,2)))
            
            c6<-removeClusters(c4,c(1,3))
            expect_equal(NCOL(clusterMatrix(c6)), 1)
            expect_equal(length(clusterType(c6)), 1)
            expect_equal(length(clusterInfo(c6)), 1)
            
            expect_error(removeClusters(c4,c(1,4)))
            c7<-removeClusters(c4,"User")
            expect_equal(NCOL(clusterMatrix(c7)), 2)
            expect_equal(length(clusterType(c7)), 2)
            expect_equal(length(clusterInfo(c7)), 2)
            
            ppC<-addClusters(cc,cbind(rep(c(-1, 1), each=5),rep(c(2, 1), each=5)),type=c("clusterMany","mergeClusters"))
            expect_equal(dim(pipelineClusters(ppC)),c(10,2))
            
            ppC<-addClusters(cc,cbind(rep(c(-1, 1), each=5)),type=c("clusterMany"))
            expect_equal(dim(pipelineClusters(ppC)),c(10,1))
            
            ppC<-addClusters(cc,cbind(rep(c(-1, 1), each=5),rep(c(2, 1), each=5)),type=c("clusterMany","mergeClusters_1"))
            expect_equal(dim(pipelineClusters(ppC)),c(10,1))
            expect_equal(dim(pipelineClusters(ppC,iteration=NA)),c(10,2))
            expect_null(pipelineClusters(cc,iteration=NA))
            
            x<-clusterLegend(cc)
            clusterLegend(cc)<-x
            clusterLegend(c4)[1]<-x
            clusterLegend(c4)[[1]]<-x[[1]]
            
            expect_error(clusterLegend(c4)[2]<-x,"must be matrix with")
            expect_error(clusterLegend(c4)[[2]]<-x[[1]],"must be matrix with")
        })
test_that("accessing transformed data works as promised", 
          {
#check all of the option handling on the dimensionality reduction arguments
  expect_equal(dim(transform(cc)), dim(assay(cc)))
  expect_equal(dim(transform(cc,dimReduce="PCA",nPCADims=3)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce="mostVar",nVarDims=3)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce=c("PCA","mostVar"),nVarDims=2)),c(2,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce=c("PCA","mostVar"),nPCADims=2)),c(2,NCOL(assay(cc))))
  expect_equal(length(transform(cc,dimReduce="mostVar",nVarDims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce="PCA",nPCADims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","mostVar"),nPCADims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","mostVar"),nVarDims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","mostVar"),nPCADims=c(2,3),nVarDims=4)),3)
  expect_equal(length(transform(cc,dimReduce=c("PCA","mostVar"),nPCADims=c(3),nVarDims=4)),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","mostVar"),nPCADims=c(2),nVarDims=c(3,4))),3)
  expect_equal(dim(transform(cc,dimReduce=c("PCA","mostVar"),nPCADims=NA,nVarDims=NA)),dim(assay(cc)))
  expect_equal(dim(transform(cc,dimReduce=c("PCA"),nPCADims=NA,nVarDims=3)),dim(assay(cc)))
  expect_equal(length(transform(cc,dimReduce=c("PCA"),nPCADims=c(NA,3),nVarDims=4)),2)
  
  
            })