context("Dimensionality Reduction")
source("create_objects.R")

test_that("makeDimReduce works as promised",{
  # check for all objects that get expected
  nDim<-3
  true3<-abs(reducedDim(sceSimDataDimRed,"PCA")[,1:nDim])

  #note: cc gives rownames to everything, so need to unname it
  expect_silent(dr3<-makeDimReduce(ceSimData,dimReduce="PCA",maxDims=nDim))
  expect_equal(unname(true3),unname(abs(reducedDim(dr3,"PCA"))))

  expect_silent(dr3<-makeDimReduce(simData,dimReduce="PCA",maxDims=nDim))
  expect_equal(true3,abs(reducedDim(dr3,"PCA")))

  expect_silent(dr3<-makeDimReduce(seSimData,dimReduce="PCA",maxDims=nDim))
  expect_equal(true3,abs(reducedDim(dr3,"PCA")))

  expect_silent(dr3<-makeDimReduce(sceSimData,dimReduce="PCA",maxDims=nDim))
  expect_equal(true3,abs(reducedDim(dr3,"PCA")))

  #check don't lose them if call on existing object
  expect_silent(dr<-makeDimReduce(sceSimDataDimRed,dimReduce="PCA",maxDims=nDim))
  expect_equal(reducedDims(sceSimDataDimRed),reducedDims(dr))

  #check with maxDims<1 (picks 4 of dimensions apparently -- never checked was correct)
  expect_silent(dr2<-makeDimReduce(cc,dimReduce="PCA",maxDims=0.5))
  expect_equal(abs(reducedDim(dr2,"PCA")[,1:4]), abs(reducedDim(dr2,"PCA")))
  
})

test_that("SingleCellFilter class works as expected",{

}
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
