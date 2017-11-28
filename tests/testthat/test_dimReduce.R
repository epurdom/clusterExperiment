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
	
	#when there are no filters
	expect_silent(scf<-SingleCellFilter(sce))
	expect_silent(filterNames(scf))
	expect_null(filterStats(scf))
	expect_error(filterStats(scf,type="Filter1"),"There are no filter statistics saved for this object")
	
	set.seed(352)
	filter<-rnorm(nrow(sce))
	filter2<-matrix(rnorm(2*nrow(sce)),ncol=2)

	#tests on adding unnamed filters
	expect_error(SingleCellFilter(sce,filterStats=filter2),"filterStats matrix must have unique column names")
	expect_error(SingleCellFilter(sce,filterStats=filter),"filterStats matrix must have unique column names")
	expect_silent(SingleCellFilter(sce,filterStats=filter,filterNames=c("Filter1")))
	expect_silent(scf1<-SingleCellFilter(sce,filterStats=filter2,filterNames=c("Filter1","Filter2")))
	expect_silent(scf2<-SingleCellFilter(sce,filterStats=filter2,filterNames=c("Filter")))
	expect_equal(scf1,scf2)
	
	#value type of replace:
	expect_silent(scf<-SingleCellFilter(sce))
	expect_silent(filterStats(scf,type="Filter1")<-filter)
	expect_equal(colnames(filterStats(scf)),"Filter1")

	expect_silent(scf<-SingleCellFilter(sce))
	expect_silent(filterStats(scf,type=c("Filter1","Filter2"))<-filter2)
	expect_equal(colnames(filterStats(scf)),c("Filter1","Filter2"))

	#matrix type of replace:	
	expect_silent(scf<-SingleCellFilter(sce))
	expect_error(filterStats(scf)<-filter,"must give matrix of values with column names")
	expect_silent(scf<-SingleCellFilter(sce))
	expect_error(filterStats(scf)<-filter2,"must give matrix of values with column names")

	####Now work with existing filters
	colnames(filter2)<-c("Filter1","Filter2")
	expect_silent(scf<-SingleCellFilter(sce,filterStats=filter2))
	expect_silent(filterNames(scf))
	expect_equal(filterStats(scf),filter2)
	expect_equal(filterStats(scf,type="Filter1"),filter2[,"Filter1"])
	expect_error(filterStats(scf,type="Myfilter"),"None of the values of 'type' argument are valid filter names")
	expect_error(filterStats(scf,type=1),"unable to find an inherited method")

	##Check Replace with existing filters in place and check actually change them, etc.
	expect_silent(scf<-SingleCellFilter(sce,filterStats=filter2))
	filter3<-filter2+1
	expect_silent(filterStats(scf) <- filter3)
	expect_equal(filterStats(scf),filter3)

	expect_silent(scf<-SingleCellFilter(sce,filterStats=filter2))
	filter4<-filter2+5
	colnames(filter4)[2]<-"Myfilter"
	expect_silent(filterStats(scf) <- filter4)
	expect_equal(filterStats(scf),cbind(filter4[,"Filter1",drop=FALSE],filter2[,"Filter2",drop=FALSE],filter4[,"Myfilter",drop=FALSE]))

	#Need more checks about replacement, etc.
	
})
test_that("filterData works as expected",{
	
	###Cutoff filter
	set.seed(352)
	filter2<-matrix(rnorm(2*nrow(sce)),ncol=2)
	colnames(filter2)<-c("Filter1","Filter2")
	expect_silent(filterStats(scf)<-filter2)
	

	tf<-filter2[,"Filter1"]>1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])

	tf<-abs(filter2[,"Filter1"])>1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1,absolute=TRUE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	tf<-abs(filter2[,"Filter1"])<1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1,keepLarge=FALSE,absolute=TRUE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	tf<-filter2[,"Filter1"]<1
	expect_silent(f1<-filterData(scf,type="Filter1",cutoff=1,keepLarge=FALSE,absolute=FALSE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	#percentile number filter
	tf<-order(filter2[,"Filter1"],decreasing=FALSE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])

	tf<-order(abs(filter2[,"Filter1"]),decreasing=FALSE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10,absolute=TRUE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	tf<-order(abs(filter2[,"Filter1"]),decreasing=TRUE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10,keepLarge=FALSE,absolute=TRUE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	tf<-order((filter2[,"Filter1"]),decreasing=TRUE)[1:10]
	expect_silent(f1<-filterData(scf,type="Filter1",percentile=10,keepLarge=FALSE,absolute=FALSE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(scf)[tf,])
	
	###Need to add test for percentile in (0,1)
	
	
})

test_that("makeFilterStats works as promised",{
	expect_silent(fs<-makeFilterStats(mat,filterStats="var"))
	expect_silent(fs<-makeFilterStats(mat,filterStats=c("mean","var")))

	#does all:
	expect_silent(fs<-makeFilterStats(mat))
	for(ii in listBuiltInFilterStats()){
		if(ii=="abscv"){
			x<-abs(sqrt(apply(mat,1,var))/apply(mat,1,mean))
			expect_equal(filterStats(fs,type=ii),unname(x))
		}
		else expect_equal(filterStats(fs,type=ii),unname(apply(mat,1,ii)))		
	}


	expect_silent(fs<-makeFilterStats(cc,filterStats="var"))
	expect_silent(fs<-makeFilterStats(cc,filterStats=c("mean","var")))
	expect_silent(fs<-makeFilterStats(cc))
	expect_silent(fs<-makeFilterStats(cc,whichClusterIgnore="Cluster1"))
	
	expect_silent(fs<-makeFilterStats(se))
	expect_silent(fs<-makeFilterStats(se,filterStats="var"))
	expect_silent(fs<-makeFilterStats(se,filterStats=c("mean","var")))
	expect_equal(sort(filterNames(fs)),sort(c("mean","var")))
	expect_silent(fs<-makeFilterStats(fs,filterStats=c("var")))
	expect_equal(sort(filterNames(fs)),sort(c("mean","var")))
	expect_silent(fs<-makeFilterStats(fs,filterStats=c("abscv")))
	expect_equal(sort(filterNames(fs)),sort(c("mean","var","abscv")))
	
	expect_silent(fs<-makeFilterStats(scf,filterStats="var"))
	expect_silent(fs<-makeFilterStats(scf,filterStats=c("mean","var")))
	expect_silent(fs<-makeFilterStats(scf))
	
	
	###Check getting out correct ones.
	contData<-simData[,1:20]
    expect_silent(temp<-makeFilterStats(contData))
	v<-apply(contData,1,var)
	expect_equal(filterStats(temp,type="var"),v)
	#should keep top 10
	expect_silent(x<-filterData(temp,percentile=10,type="var",keepLarge=TRUE,absolute=FALSE))
	whTop<-order(v,decreasing=TRUE)[1:10]
	expect_equal(v[whTop],filterStats(x,type="var"))
	expect_equal(contData[whTop,],unname(assay(x)))

	###Check getting out correct ones.
	countData<-simCount[,1:20]
    expect_silent(temp<-makeFilterStats(countData,isCount=TRUE))
	v<-apply(log2(countData+1),1,var)
	expect_equal(filterStats(temp,type="var"),v)
	#should keep top 10
	expect_silent(x<-filterData(temp,percentile=10,type="var",keepLarge=TRUE,absolute=FALSE))
	whTop<-order(v,decreasing=TRUE)[1:10]
	expect_equal(v[whTop],filterStats(x,type="var"))
	expect_equal(log2(countData[whTop,]+1),unname(transformData(x,isCount=TRUE)))
	
	
})
  #
  # expect_equal(dim(transformData(cc,dimReduce="var",nFilter=3)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="abscv",nFilter=3)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="mad",nFilter=3)), c(3,NCOL(assay(cc))))
  #
  # expect_equal(dim(transformData(cc,dimReduce="PCA",nPCADims=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="var",nFilter=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="abscv",nFilter=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce="mad",nFilter=3,ignoreUnassigned=TRUE)), c(3,NCOL(assay(cc))))
  #
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var"),nFilter=2)),c(2,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var"),nPCADims=4)),c(4,NCOL(assay(cc))))
  # expect_equal(length(transformData(cc,dimReduce="var",nFilter=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce="PCA",nPCADims=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nFilter=c(2,3))),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(2,3),nFilter=4)),3)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(3),nFilter=4)),2)
  # expect_equal(length(transformData(cc,dimReduce=c("PCA","var"),nPCADims=c(2),nFilter=c(3,4))),3)
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var"),nPCADims=NA,nFilter=NA)),dim(assay(cc)))
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA"),nPCADims=NA,nFilter=3)),dim(assay(cc)))
  # expect_equal(length(transformData(cc,dimReduce=c("PCA"),nPCADims=c(NA,3),nFilter=4)),2)
  #
  # expect_equal(length(transformData(cc,dimReduce=c("var","abscv","mad"),nPCADims=c(NA,3),nFilter=4)),3)
  # expect_equal(length(transformData(cc,dimReduce=c("var","abscv","mad"),nPCADims=c(NA,3),nFilter=c(2,4))),6)
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA","var","abscv"),nPCADims=c(3),nFilter=NA)),c(3,NCOL(assay(cc))))
  # expect_equal(dim(transformData(cc,dimReduce=c("PCA"),nPCADims=c(3),nFilter=2)),c(3,NCOL(assay(cc))))
