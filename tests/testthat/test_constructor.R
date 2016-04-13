library(clusterExperiment)

mat <- matrix(data=rnorm(200), ncol=10)
mat[1,1]<- -1 #force a negative value
labels <- as.character(gl(5, 2))
labels[c(1:2)]<- c("-1","-2") #make sure some not assigned
labels<-factor(labels)
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
                         "must be a vector of length equal")
            expect_error(clusterExperiment(as.data.frame(mat), labels, function(x){x}),
                         "unable to find an inherited method for function")
            expect_warning(clusterExperiment(mat, as.character(labels), function(x){x}),"was coerced to integer values")
            expect_warning(clusterExperiment(mat, labels, function(x){x}), "was coerced to integer values")

            expect_is(cc, "ClusterExperiment")
            expect_is(cc, "SummarizedExperiment")
            
            expect_equal(nSamples(cc),ncol(mat))
            expect_equal(nFeatures(cc),nrow(mat))
            expect_equal(nClusters(cc),1)

            clusterExperiment(se,labMat,transformation=function(x){x})
            expect_warning(clusterExperiment(se,labels,transformation=function(x){x}))
            expect_warning(clusterExperiment(se,as.character(labels),transformation=function(x){x}))
            
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
 
            expect_null(clusterLabels(c3))
            clusterLabels(c3)<-c("User","User")
            expect_equal(length(clusterLabels(c3)),2)
            expect_error(clusterLabels(c3,"pipeline"))
            
              #check adding matrix of clusters
            c4<-addClusters(cc,cbind(rep(c(-1, 1), each=5),rep(c(2, 1), each=5)),type="New")
            expect_equal(NCOL(clusterMatrix(c4)), 3)
            expect_equal(length(clusterType(c4)), 3)
            expect_equal(length(clusterInfo(c4)), 3)
            expect_equal(primaryCluster(c4), primaryCluster(cc))
            
            ccR<-cc
            dendrogram(ccR)<-NULL
            expect_equal(ccR, removeUnclustered(cc))

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
            
            x<-clusterColors(cc)
            clusterColors(cc)<-x
            clusterColors(c4)[1]<-x
            clusterColors(c4)[[1]]<-x[[1]]
            
            expect_error(clusterColors(c4)[2]<-x,"must be matrix with")
            expect_error(clusterColors(c4)[[2]]<-x[[1]],"must be matrix with")
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