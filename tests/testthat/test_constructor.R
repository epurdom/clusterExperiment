context("Constructor")
source("create_objects.R")

test_that("`clusterExperiment` constructor works with matrix and
          SummarizedExperiments", {
            expect_error(clusterExperiment(mat), "missing")
            expect_error(clusterExperiment(mat,as.numeric(numLabels),transformation=log),info="Error checking transFun")
            print(numLabels)
            expect_error(clusterExperiment(mat, as.numeric(numLabels)), "missing")
            expect_error(clusterExperiment(mat, numLabels[1:2], function(x){x}),
                         "must be a matrix of rows equal")
            expect_error(clusterExperiment(as.data.frame(mat), numLabels, function(x){x}),
                         "unable to find an inherited method for function")
            #test character input
            ccChar<-clusterExperiment(mat, chLabels, function(x){x})
            expect_is(primaryCluster(ccChar),"numeric")
            expect_is(primaryClusterNamed(ccChar),"character")
            expect_equal(sort(unique(primaryClusterNamed(ccChar))),sort(unique(chLabels)))

            #test factor input
            clusterExperiment(mat, numLabels, function(x){x})

            expect_is(cc, "ClusterExperiment")
            expect_is(cc, "SummarizedExperiment")

            expect_equal(nSamples(cc),ncol(mat))
            expect_equal(nFeatures(cc),nrow(mat))
            expect_equal(nClusters(cc),2)
            expect_equal(NCOL(clusterMatrix(ccSE)), NCOL(labMat))
            
            expect_equal(colData(ccSE),colData(se)) 
            expect_equal(rownames(ccSE),rownames(se)) 
            expect_equal(colnames(ccSE),colnames(se)) 
            expect_equal(metadata(ccSE),metadata(se)) 
            expect_equal(rowData(ccSE),rowData(se)) 
            
            show(ccSE)
 })

test_that("adding clusters, setting primary labels and remove unclustered cells
          work as promised", {
            
            expect_is(transformation(ccSE),"function")

            ##########
            #subsetting 
            #right now just checks SE info (which isn't biggest place needs unit tests since uses callNextMethod() so should work)
            #therefore currently just really checks no errors. 
            x1<-ccSE[,5:8]
            ###Check retain SE info
            expect_equal(colData(x1),colData(se[,5:8]) )
            expect_equal(rownames(x1),rownames(se[,5:8])) 
            expect_equal(colnames(x1),colnames(se[,5:8])) 
            expect_equal(metadata(x1),metadata(se[,5:8])) 
            expect_equal(rowData(x1),rowData(se[,5:8])) 
            
            x2<-ccSE[1:2,]
            ###Check retain SE info
            expect_equal(colData(x2),colData(se[1:2,]) )
            expect_equal(rownames(x2),rownames(se[1:2,])) 
            expect_equal(colnames(x2),colnames(se[1:2,])) 
            expect_equal(metadata(x2),metadata(se[1:2,])) 
            expect_equal(rowData(x2),rowData(se[1:2,])) 
            
            x3<-ccSE[,3]
            ###Check retain SE info
            expect_equal(colData(x3),colData(se[,3]) )
            expect_equal(rownames(x3),rownames(se[,3])) 
            expect_equal(colnames(x3),colnames(se[,3])) 
            expect_equal(metadata(x3),metadata(se[,3])) 
            expect_equal(rowData(x3),rowData(se[,3])) 
            
            x4<-ccSE[3,]
            ###Check retain SE info
            expect_equal(colData(x4),colData(se[3,]) )
            expect_equal(rownames(x4),rownames(se[3,])) 
            expect_equal(colnames(x4),colnames(se[3,])) 
            expect_equal(metadata(x4),metadata(se[3,])) 
            expect_equal(rowData(x4),rowData(se[3,])) 
            
            ##########
            #addClusters
            c1 <- addClusters(ccSE, rep(c(-1, 1, 2), each=5),clusterTypes="newUser")
            ###Check retain SE info
            expect_equal(colData(c1),colData(se)) 
            expect_equal(rownames(c1),rownames(se)) 
            expect_equal(colnames(c1),colnames(se)) 
            expect_equal(metadata(c1),metadata(se)) 
            expect_equal(rowData(c1),rowData(se)) 
            #Other checks
            expect_equal(NCOL(clusterMatrix(c1)), nClusters(ccSE)+1)
            expect_equal(unname(clusterTypes(c1)), unname(c(clusterTypes(ccSE),"newUser")))
            expect_equal(length(clusterInfo(c1)), nClusters(ccSE)+1)
            expect_equal(primaryCluster(c1), primaryCluster(ccSE))
            primaryClusterIndex(c1) <- 3
            expect_false(all(primaryCluster(c1)==primaryCluster(ccSE)))

            ####check adding a clusterExperiment to existing CE
            expect_error(addClusters(ccSE,smSimCE),"Cannot merge clusters from different data") #assays don't match
            c3<-addClusters(ccSE,ccSE)
            expect_equal(NCOL(clusterMatrix(c3)), nClusters(ccSE)*2)
            expect_equal(length(clusterTypes(c3)), nClusters(ccSE)*2)
            expect_equal(length(clusterInfo(c3)), nClusters(ccSE)*2)
            expect_equal(primaryCluster(c3), primaryCluster(ccSE))
            ###Check retain SE info
            expect_equal(colData(c3),colData(se)) 
            expect_equal(rownames(c3),rownames(se)) 
            expect_equal(colnames(c3),colnames(se)) 
            expect_equal(metadata(c3),metadata(se)) 
            expect_equal(rowData(c3),rowData(se)) 
            #Other checks after adding CE object  
            expect_error(clusterLabels(c3)[1:2]<-c("User","User"),"cannot have duplicated clusterLabels")
            clusterLabels(c3)[1:2]<-c("User1","User2")
            clusterLabels(c3)[1]<-"User4"
            expect_error(clusterLabels(c3)[1]<-"User2","duplicated clusterLabels")
            expect_equal(length(clusterLabels(c3)),nClusters(ccSE)*2)
            
            ###check adding matrix of clusters
            c4<-addClusters(ccSE,clusterMatrix(ccSE),clusterTypes="New")
            newLeng<-2*nClusters(ccSE)
            expect_equal(NCOL(clusterMatrix(c4)), newLeng)
            expect_equal(length(clusterTypes(c4)), newLeng)
            expect_equal(length(clusterInfo(c4)), newLeng)
            expect_equal(primaryCluster(c4), primaryCluster(ccSE))
            ###Check retain SE info
            expect_equal(colData(c4),colData(se)) 
            expect_equal(rownames(c4),rownames(se)) 
            expect_equal(colnames(c4),colnames(se)) 
            expect_equal(metadata(c4),metadata(se)) 
            expect_equal(rowData(c4),rowData(se)) 
    
            
            ##########
            #check removing unclustered
            
            #no -1 in primary cluster
            matTemp<-abs(labMat)
            ccTemp<-clusterExperiment(mat,matTemp,transformation=function(x){x})
            expect_equal(ccTemp, removeUnclustered(ccTemp)) 
            
            #-1 in primary cluster
            whUn<-which(primaryCluster(ccSE) <0)
            ccR<-removeUnclustered(ccSE)
            expect_equal(NCOL(ccR), NCOL(ccSE)-length(whUn))
            ###Check retain SE info
            expect_equal(colData(ccR),colData(se[,-whUn]) )
            expect_equal(rownames(ccR),rownames(se)) 
            expect_equal(colnames(ccR),colnames(se[,-whUn])) 
            expect_equal(metadata(ccR),metadata(se)) 
            expect_equal(rowData(ccR),rowData(se)) 

            ##########
            #check removeClusters
            
            #single cluster
            c5<-removeClusters(c4,1)
            expect_equal(NCOL(clusterMatrix(c5)), nClusters(c4)-1)
            expect_equal(length(clusterTypes(c5)), nClusters(c4)-1)
            expect_equal(length(clusterInfo(c5)), nClusters(c4)-1)
            expect_equal(primaryCluster(c4), primaryCluster(removeClusters(c4,2)))
            ###Check retain SE info 
            expect_equal(colData(c5),colData(se)) 
            expect_equal(rownames(c5),rownames(se)) 
            expect_equal(colnames(c5),colnames(se)) 
            expect_equal(metadata(c5),metadata(se)) 
            expect_equal(rowData(c5),rowData(se)) 
            
            #vector clusters
            c6<-removeClusters(c4,c(1,3))
            expect_equal(NCOL(clusterMatrix(c6)), nClusters(c4)-2)
            expect_equal(length(clusterTypes(c6)), nClusters(c4)-2)
            expect_equal(length(clusterInfo(c6)), nClusters(c4)-2)
            ###Check retain SE info 
            expect_equal(colData(c6),colData(se)) 
            expect_equal(rownames(c6),rownames(se)) 
            expect_equal(colnames(c6),colnames(se)) 
            expect_equal(metadata(c6),metadata(se)) 
            expect_equal(rowData(c6),rowData(se)) 
            
            expect_error(removeClusters(c4,c(1,nClusters(c4)+1)),"invalid indices")
            
            c7<-removeClusters(c4,"User") #two have "user" label
            expect_equal(NCOL(clusterMatrix(c7)), nClusters(c4)-2)
            expect_equal(length(clusterTypes(c7)), nClusters(c4)-2)
            expect_equal(length(clusterInfo(c7)), nClusters(c4)-2)

            ##########
            #check workflow stuff
            ppC<-addClusters(cc,cbind(rep(c(-1, 1,2), each=5),rep(c(2, 1,3), each=5)),clusterTypes=c("clusterMany","mergeClusters"))
            expect_equal(dim(workflowClusters(ppC)),c(nSamples(cc),2))

            ppC<-addClusters(cc,cbind(rep(c(-1, 1,2), each=5)),clusterTypes=c("clusterMany"))
            expect_equal(dim(workflowClusters(ppC)),c(nSamples(cc),1))

            ppC<-addClusters(cc,cbind(rep(c(-1, 1,2), each=5),rep(c(2, 3,1), each=5)),clusterTypes=c("clusterMany","mergeClusters.1"))
            expect_equal(dim(workflowClusters(ppC)),c(nSamples(cc),1))
            expect_equal(dim(workflowClusters(ppC,iteration=NA)),c(nSamples(cc),2))
            expect_null(workflowClusters(cc,iteration=NA))

            ##########
            #clusterLegend
            x<-clusterLegend(cc)
            clusterLegend(cc)<-x
            clusterLegend(c4)[1:2]<-x[1:2]
            clusterLegend(c4)[[1]]<-x[[1]]
#add wrong dimensions:
            expect_error(clusterLegend(c4)[3]<-list(x[[1]][1:2,]),"each element of `clusterLegend` must be matrix with")
            expect_error(clusterLegend(c4)[[3]]<-x[[1]][1:2,],"must be matrix with")
            
        })
test_that("accessing transformed data works as promised",
          {
#check all of the option handling on the dimensionality reduction arguments
  expect_equal(dim(transform(cc)), dim(assay(cc)))
  expect_equal(dim(transform(cc,dimReduce="PCA",nPCADims=3)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce="PCA",nPCADims=0.5)), c(4,NCOL(assay(cc))))

  expect_equal(dim(transform(cc,dimReduce="PCA",nPCADims=c(8,0.5,3))[[2]]), c(4,NCOL(assay(cc))))

  expect_equal(dim(transform(cc,dimReduce="var",nVarDims=3)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce="cv",nVarDims=3)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce="mad",nVarDims=3)), c(3,NCOL(assay(cc))))

  expect_equal(dim(transform(cc,dimReduce="PCA",nPCADims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce="var",nVarDims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce="cv",nVarDims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce="mad",nVarDims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  
  expect_equal(dim(transform(cc,dimReduce=c("PCA","var"),nVarDims=2)),c(2,NCOL(assay(cc))))
  expect_equal(dim(transform(cc,dimReduce=c("PCA","var"),nPCADims=4)),c(4,NCOL(assay(cc))))
  expect_equal(length(transform(cc,dimReduce="var",nVarDims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce="PCA",nPCADims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","var"),nPCADims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","var"),nVarDims=c(2,3))),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","var"),nPCADims=c(2,3),nVarDims=4)),3)
  expect_equal(length(transform(cc,dimReduce=c("PCA","var"),nPCADims=c(3),nVarDims=4)),2)
  expect_equal(length(transform(cc,dimReduce=c("PCA","var"),nPCADims=c(2),nVarDims=c(3,4))),3)
  expect_equal(dim(transform(cc,dimReduce=c("PCA","var"),nPCADims=NA,nVarDims=NA)),dim(assay(cc)))
  expect_equal(dim(transform(cc,dimReduce=c("PCA"),nPCADims=NA,nVarDims=3)),dim(assay(cc)))
  expect_equal(length(transform(cc,dimReduce=c("PCA"),nPCADims=c(NA,3),nVarDims=4)),2)

  expect_equal(length(transform(cc,dimReduce=c("var","cv","mad"),nPCADims=c(NA,3),nVarDims=4)),3)
  expect_equal(length(transform(cc,dimReduce=c("var","cv","mad"),nPCADims=c(NA,3),nVarDims=c(2,4))),6)  
  expect_equal(dim(transform(cc,dimReduce=c("PCA","var","cv"),nPCADims=c(3),nVarDims=NA)),c(3,NCOL(assay(cc))))  
  expect_equal(dim(transform(cc,dimReduce=c("PCA"),nPCADims=c(3),nVarDims=2)),c(3,NCOL(assay(cc))))  
          })

test_that("workflow functions work",
          {
  ceNew<-combineMany(ceSim,proportion=0.7)
  ceNew<-combineMany(ceNew,proportion=0.3,clusterLabel="combineMany,v2")
  expect_equal(clusterLabels(ceNew)[1:2],c("combineMany,v2","combineMany.1"))
  expect_equal(clusterTypes(ceNew)[1:2],c("combineMany","combineMany.1"))
  ceNew2<-setToCurrent(ceNew,whichCluster="combineMany.1")
  expect_equal(clusterLabels(ceNew2)[1:2],c("combineMany,v2","combineMany"))
  expect_equal(clusterTypes(ceNew2)[1:2],c("combineMany.2","combineMany"))
  ceNew3<-setToCurrent(ceNew2,whichCluster="combineMany.2")
  expect_equal(clusterLabels(ceNew3)[1:2],c("combineMany,v2","combineMany.3"))
  expect_equal(clusterTypes(ceNew3)[1:2],c("combineMany","combineMany.3"))

  ceNew4<-setToFinal(ceNew,whichCluster="combineMany,v2",clusterLabel="Final Version")
  expect_equal(primaryClusterIndex(ceNew4),1)
  expect_equal(clusterLabels(ceNew4)[primaryClusterIndex(ceNew4)],"Final Version")
  expect_equal(clusterTypes(ceNew4)[primaryClusterIndex(ceNew4)],"final")
  
  ceNew5<-setToFinal(ceNew,whichCluster="combineMany.1",clusterLabel="Final Version")
  expect_equal(primaryClusterIndex(ceNew5),2)
  expect_equal(clusterLabels(ceNew5)[primaryClusterIndex(ceNew5)],"Final Version")
  expect_equal(clusterTypes(ceNew5)[primaryClusterIndex(ceNew5)],"final")
          })