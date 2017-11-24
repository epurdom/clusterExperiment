context("Dimensionality Reduction")
source("create_objects.R")

test_that("makeDimReduce works as promised",{

  #check makeDimReduce with CE object
  expect_silent(dr<-makeDimReduce(cc,dimReduce="PCA",maxDims=3))
  expect_equal(dim(reducedDim(dr,"PCA")), c(NCOL(assay(cc)),3))
  expect_silent(dr2<-makeDimReduce(cc,dimReduce="PCA",maxDims=0.5))
  expect_equal(dim(reducedDim(dr2,"PCA")), c(NCOL(assay(cc)),4))
  
  # check for all objects that get expected
  true3<-reducedDim(sceSimDataDimRed,"PCA")[,1:3]
  expect_silent(dr3<-makeDimReduce(simData,dimReduce="PCA",maxDims=3))
  expect_equal(true3,reducedDim(dr3,"PCA"))
  expect_silent(dr3<-makeDimReduce(seSimData,dimReduce="PCA",maxDims=3))
  expect_equal(true3,reducedDim(dr3,"PCA"))
  expect_silent(dr3<-makeDimReduce(sceSimData,dimReduce="PCA",maxDims=3))
  expect_equal(true3,reducedDim(dr3,"PCA"))

  #check don't lose them
  expect_silent(dr<-makeDimReduce(sceSimDataDimRed,dimReduce="PCA",maxDims=3))
  expect_equal(reducedDim(sceSimDataDimRed,"PCA"),reducedDim(dr,"PCA"))
  
})
  # expect_equal(dim(transformData(cc,dimReduce="PCA",nPCADims=c(8,0.5,3))[[2]]), c(4,NCOL(assay(cc))))
  #
  # expect_equal(dim(transformData(cc,dimReduce="var",nVarDims=3)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="cv",nVarDims=3)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="mad",nVarDims=3)), c(3,NCOL(assay(cc))))
  #
  # expect_equal(dim(transformData(cc,dimReduce="PCA",nPCADims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="var",nVarDims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="cv",nVarDims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="mad",nVarDims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  #
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var"),nVarDims=2)),c(2,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var"),nPCADims=4)),c(4,NCOL(assay(cc))))
  # expect_equal(length(transformData(cc,dimReduce="var",nVarDims=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce="PCA",nPCADims=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nVarDims=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(2,3),nVarDims=4)),3)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(3),nVarDims=4)),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(2),nVarDims=c(3,4))),3)
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var"),nPCADims=NA,nVarDims=NA)),dim(assay(cc)))
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA"),nPCADims=NA,nVarDims=3)),dim(assay(cc)))
  # expect_equal(length(transformData(cc,dimReduce=c("PCA"),nPCADims=c(NA,3),nVarDims=4)),2)
  #
  # expect_equal(length(transformData(cc,dimReduce=c("var","cv","mad"),nPCADims=c(NA,3),nVarDims=4)),3)
  # expect_equal(length(transformData(cc,dimReduce=c("var","cv","mad"),nPCADims=c(NA,3),nVarDims=c(2,4))),6)
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var","cv"),nPCADims=c(3),nVarDims=NA)),c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA"),nPCADims=c(3),nVarDims=2)),c(3,NCOL(assay(cc))))
