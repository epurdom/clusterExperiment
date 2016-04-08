library(clusterCells)

mat <- matrix(data=rnorm(200), ncol=10)
mat[1,1]<- -1 #force a negative value
labels <- gl(5, 2)
se <- SummarizedExperiment(mat)

cc <- clusterCells(mat, as.numeric(labels), transformation = function(x){x})
cc2 <- clusterCells(se, as.numeric(labels), transformation = function(x){x})

test_that("`clusterCells` constructor works with matrix and
          SummarizedExperiments", {
            expect_equal(cc, cc2)
            expect_error(clusterCells(mat), "missing")
            expect_error(clusterCells(mat,as.numeric(labels),transformation=log),info="Error checking transFun")
            expect_error(clusterCells(mat, as.numeric(labels)), "missing")
            expect_error(clusterCells(mat, labels[1:2], function(x){x}),
                         "must be a vector of length equal to the number of samples")
            expect_error(clusterCells(as.data.frame(mat), labels, function(x){x}),
                         "must be a matrix or SummarizedExperiment object")
            expect_error(clusterCells(mat, as.character(labels), function(x){x}),
                         "must be a numeric matrix")

            expect_warning(clusterCells(mat, labels, function(x){x}), "was coerced to numeric")

            expect_is(cc, "ClusterCells")
            expect_is(cc, "SummarizedExperiment")
          })

test_that("adding clusters, setting primary labels and remove unclustered cells
          work as promised", {
            expect_equal(NCOL(allClusters(cc)), 1)
            expect_is(transformation(cc),"function")

            c1 <- addClusters(cc, rep(c(-1, 1), each=5))
            expect_equal(NCOL(allClusters(c1)), 2)
            expect_equal(length(clusterType(c1)), 2)
            expect_equal(length(clusterInfo(c1)), 2)
            expect_equal(primaryCluster(c1), primaryCluster(cc))
            primaryClusterIndex(c1) <- 2
            expect_false(all(primaryCluster(c1)==primaryCluster(cc)))

              #check adding a clusterCellsObject
            c3<-addClusters(cc,cc)
            expect_equal(NCOL(allClusters(c3)), 2)
            expect_equal(length(clusterType(c3)), 2)
            expect_equal(length(clusterInfo(c3)), 2)
            expect_equal(primaryCluster(c3), primaryCluster(cc))
 
              #check adding matrix of clusters
            c4<-addClusters(cc,cbind(rep(c(-1, 1), each=5),rep(c(2, 1), each=5)),type="New")
            expect_equal(NCOL(allClusters(c4)), 3)
            expect_equal(length(clusterType(c4)), 3)
            expect_equal(length(clusterInfo(c4)), 3)
            expect_equal(primaryCluster(c4), primaryCluster(cc))
            
            expect_equal(cc, removeUnclustered(cc))

            c2 <- removeUnclustered(c1)
            expect_equal(NCOL(c2), 5)
            
            #check removing index of clusters
            c5<-removeClusters(c4,1)
            expect_equal(NCOL(allClusters(c5)), 2)
            expect_equal(length(clusterType(c5)), 2)
            expect_equal(length(clusterInfo(c5)), 2)
            expect_equal(primaryCluster(c4), primaryCluster(removeClusters(c4,2)))
            
            c6<-removeClusters(c4,c(1,3))
            expect_equal(NCOL(allClusters(c6)), 1)
            expect_equal(length(clusterType(c6)), 1)
            expect_equal(length(clusterInfo(c6)), 1)
            
            expect_error(removeClusters(c4,c(1,4)))
            c7<-removeClusters(c4,"User")
            expect_equal(NCOL(allClusters(c7)), 2)
            expect_equal(length(clusterType(c7)), 2)
            expect_equal(length(clusterInfo(c7)), 2)
            
            ppC<-addClusters(cc,cbind(rep(c(-1, 1), each=5),rep(c(2, 1), each=5)),type=c("compareChoices","mergeClusters"))
            expect_equal(dim(pipelineClusters(ppC)),c(10,2))
            
            ppC<-addClusters(cc,cbind(rep(c(-1, 1), each=5)),type=c("compareChoices"))
            expect_equal(dim(pipelineClusters(ppC)),c(10,1))
            
            ppC<-addClusters(cc,cbind(rep(c(-1, 1), each=5),rep(c(2, 1), each=5)),type=c("compareChoices","mergeClusters_1"))
            expect_equal(dim(pipelineClusters(ppC)),c(10,1))
            expect_equal(dim(pipelineClusters(ppC,iteration=NA)),c(10,2))
            expect_null(pipelineClusters(cc,iteration=NA))
            
            
        })
test_that("accessing transformed data works as promised", {
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
  
            })