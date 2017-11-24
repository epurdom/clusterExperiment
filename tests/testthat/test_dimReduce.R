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
	expect_silent(scf<-SingleCellFilter(sce))
	expect_silent(filterNames(scf))
	expect_null(filterStats(scf))
	expect_error(filterStat(scf,type="Filter1"),"There are no filter statistics saved for this object")
	
	filter<-rnorm(nrow(sce))
	expect_silent(SingleCellFilter(sce,filterStats=filter))
	
	set.seed(352)
	filter<-matrix(rnorm(2*nrow(sce)),ncol=2)
	expect_error(SingleCellFilter(sce,filterStats=filter),"filterStats matrix must have unique column names")
	colnames(filter)<-c("Filter1","Filter2")
	expect_silent(scf<-SingleCellFilter(sce,filterStats=filter))
	
	expect_silent(filterNames(scf))
	expect_equal(filterStats(scf),filter)
	expect_equal(filterStat(scf,type="Filter1"),filter[,"Filter1"])
	expect_error(filterStat(scf,type="Myfilter"),"is not the name of a filter statistic held by the object")
	expect_error(filterStat(scf,type=1),"unable to find an inherited method")
	
	###Cutoff filter
	tf<-filter[,"Filter1"]>1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])

	tf<-abs(filter[,"Filter1"])>1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1,absolute=TRUE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	tf<-abs(filter[,"Filter1"])<1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1,keepLarge=FALSE,absolute=TRUE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	tf<-filter[,"Filter1"]<1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1,keepLarge=FALSE,absolute=FALSE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	#percentile number filter
	tf<-order(filter[,"Filter1"],decreasing=FALSE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])

	tf<-order(abs(filter[,"Filter1"]),decreasing=FALSE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10,absolute=TRUE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	tf<-order(abs(filter[,"Filter1"]),decreasing=TRUE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10,keepLarge=FALSE,absolute=TRUE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	tf<-order((filter[,"Filter1"]),decreasing=TRUE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10,keepLarge=FALSE,absolute=FALSE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	###Need to add test for percentile in (0,1)
	
	
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
