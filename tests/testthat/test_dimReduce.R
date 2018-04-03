context("Dimensionality Reduction")


test_that("makeReducedDims works as promised",{
  # check for all objects that get expected
  nDim<-3
  true3<-abs(reducedDim(sceSimDataDimRed,"PCA")[,1:nDim])

  #note: cc gives rownames to everything, so need to unname it
  expect_silent(dr3<-makeReducedDims(ceSimData,reducedDims="PCA",maxDims=nDim))
  expect_equal(unname(true3),unname(abs(reducedDim(dr3,"PCA"))))

  expect_silent(dr3<-makeReducedDims(simData,reducedDims="PCA",maxDims=nDim))
  expect_equal(true3,abs(reducedDim(dr3,"PCA")))

  expect_silent(dr3<-makeReducedDims(seSimData,reducedDims="PCA",maxDims=nDim))
  expect_equal(true3,abs(reducedDim(dr3,"PCA")))

  expect_silent(dr3<-makeReducedDims(sceSimData,reducedDims="PCA",maxDims=nDim))
  expect_equal(true3,abs(reducedDim(dr3,"PCA")))

  #check don't lose them if call on existing object
  expect_silent(dr<-makeReducedDims(sceSimDataDimRed,reducedDims="PCA",maxDims=nDim))
  expect_equal(reducedDims(sceSimDataDimRed),reducedDims(dr))

  #check with maxDims<1 (picks 4 of dimensions apparently -- never checked was correct)
  expect_silent(dr2<-makeReducedDims(cc,reducedDims="PCA",maxDims=0.5))
  expect_equal(abs(reducedDim(dr2,"PCA")[,1:4]), abs(reducedDim(dr2,"PCA")))
  
})

test_that("Filter functions work as expected",{
	
	#when there are no filters
	expect_silent(filterNames(sceSimData))
	expect_equal(clusterExperiment:::filterStats(sceSimData),rowData(sceSimData)[,"b",drop=FALSE])
	expect_error(clusterExperiment:::filterStats(sceSimData,type="Filter1"),"None of the values of 'type' argument are valid filter names")
	
	##Redo when no rowData.
	testSe<-SummarizedExperiment(mat)
	expect_equal(length(filterNames(testSe)),0)
	expect_equal(ncol(clusterExperiment:::filterStats(testSe)),0)
	expect_error(clusterExperiment:::filterStats(testSe,type="Filter1"),"There are no filter statistics saved for this object")
	
	##Redo when no valid filters.
	testSe<-SummarizedExperiment(mat,rowData=gData[,c("a","c")])
	expect_equal(length(filterNames(testSe)),0)
	expect_equal(ncol(clusterExperiment:::filterStats(testSe)),0)
	expect_error(clusterExperiment:::filterStats(testSe,type="Filter1"),"None of the values of 'type' argument are valid filter names")
	
	set.seed(352)
	filter<-rnorm(nrow(sce))
	filter2<-matrix(rnorm(2*nrow(sce)),ncol=2)

	#tests on adding unnamed filters	
	#value type of replace:
	sceTest<-sce
	expect_silent(clusterExperiment:::filterStats(sceTest,type="Filter1")<-filter)
	expect_equal(colnames(clusterExperiment:::filterStats(sceTest)),c("b","Filter1"))

	sceTest<-sce
	expect_silent(clusterExperiment:::filterStats(sceTest,type=c("Filter1","Filter2"))<-filter2)
	expect_equal(colnames(clusterExperiment:::filterStats(sceTest)),c("b","Filter1","Filter2"))

	#matrix type of replace:	
	sceTest<-sce
	expect_error(clusterExperiment:::filterStats(sceTest)<-filter,"must give matrix of values with column names")
	sceTest<-sce
	expect_error(clusterExperiment:::filterStats(sceTest)<-filter2,"must give matrix of values with column names")

	####Now work with existing filters
	colnames(filter2)<-c("Filter1","Filter2")
	sceTest<-sce
	expect_silent(clusterExperiment:::filterStats(sceTest)<-filter2)
	expect_silent(filterNames(sceTest))
	expect_equal(clusterExperiment:::filterStats(sceTest),DataFrame(b=rowData(sce)[,"b"],filter2))
	expect_equal(clusterExperiment:::filterStats(sceTest,type="Filter1"),DataFrame(filter2)[,"Filter1",drop=FALSE])
	expect_error(clusterExperiment:::filterStats(sceTest,type="Myfilter"),"None of the values of 'type' argument are valid filter names")
	expect_error(clusterExperiment:::filterStats(sceTest,type=1),"unable to find an inherited method")

	##Check Replace with existing filters in place and check actually change them, etc.
	sceTest<-SummarizedExperiment(assay(sce))
	expect_silent(clusterExperiment:::filterStats(sceTest)<-filter2)
	filter3<-filter2+1
	expect_silent(clusterExperiment:::filterStats(sceTest) <- filter3)
	expect_equal(clusterExperiment:::filterStats(sceTest),DataFrame(filter3))

	sceTest<-SummarizedExperiment(assay(sce))
	expect_silent(clusterExperiment:::filterStats(sceTest)<-filter2)
	filter4<-filter2+5
	colnames(filter4)[2]<-"Myfilter"
	expect_silent(clusterExperiment:::filterStats(sceTest) <- filter4)
	expect_equal(clusterExperiment:::filterStats(sceTest),DataFrame(filter4[,"Filter1",drop=FALSE],filter2[,"Filter2",drop=FALSE],filter4[,"Myfilter",drop=FALSE]))

	#Need more checks about replacement, etc.
	
})

test_that("reduce and filter work with hdf5",{
	expect_silent(filterNames(hdfSCE))

	expect_silent(f1<-filterData(hdfSCE,filterStats="Filter1",cutoff=1))
	expect_silent(f1<-filterData(hdfSCE,filterStats="Filter1",cutoff=1))
	
	expect_silent(fs<-makeFilterStats(hdfSCE,filterStats="var"))
	expect_silent(fs<-makeFilterStats(hdfSCE,filterStats=c("mean","var")))
	expect_silent(out<-filterData(fs,filterStats=c("mean"),cutoff=1))


	expect_silent(defaultNDims(hdfObj,"PCA"))

	#add pca to it
	nDim<-3
	expect_silent(dr3<-makeReducedDims(hdfObj,reducedDims="PCA",maxDims=nDim))
	expect_equal(defaultNDims(dr3,"PCA"),3)

	#test transformation -- need make CE object
    expect_silent(clustNothing1 <- clusterSingle(hdfObj,
		  reduceMethod = "none",
  		  mainClusterArgs=list( clusterArgs=list(k=3),clusterFunction=listBuiltInFunctions()[[1]]),
  	       subsample=FALSE, sequential=FALSE, isCount=FALSE))  	
	
	transformation(clustNothing1)<-function(x){exp(x)}
	expect_equal(exp(assay(clustNothing1)),unname(transformData(clustNothing1)))
	

	
})


test_that("filterData works as expected",{
	
	###Cutoff filter
	set.seed(352)
	filter2<-matrix(rnorm(2*nrow(sceSimData)),ncol=2)
	colnames(filter2)<-c("Filter1","Filter2")
	expect_silent(clusterExperiment:::filterStats(sceSimData)<-filter2)
	

	tf<-filter2[,"Filter1"]>1
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",cutoff=1))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])

	tf<-abs(filter2[,"Filter1"])>1
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",cutoff=1,absolute=TRUE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])
	
	tf<-abs(filter2[,"Filter1"])<1
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",cutoff=1,keepLarge=FALSE,absolute=TRUE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])
	tf<-filter2[,"Filter1"]<1
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",cutoff=1,keepLarge=FALSE,absolute=FALSE))
	expect_equal(NROW(f1),sum(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])
	
	#percentile number filter
	tf<-order(filter2[,"Filter1"],decreasing=TRUE)[1:10]
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",percentile=10))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])

	tf<-order(abs(filter2[,"Filter1"]),decreasing=TRUE)[1:10]
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",percentile=10,absolute=TRUE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])
	
	tf<-order(abs(filter2[,"Filter1"]),decreasing=FALSE)[1:10]
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",percentile=10,keepLarge=FALSE,absolute=TRUE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])
	tf<-order((filter2[,"Filter1"]),decreasing=FALSE)[1:10]
	expect_silent(f1<-filterData(sceSimData,filterStats="Filter1",percentile=10,keepLarge=FALSE,absolute=FALSE))
	expect_equal(NROW(f1),length(tf))
	expect_equal(assay(f1),assay(sceSimData)[tf,])
	
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
			x<-unname(x)
		}
		else x<-(unname(apply(mat,1,ii)))
		x<-DataFrame(x)
		colnames(x)<-ii
		expect_equal(clusterExperiment:::filterStats(fs,type=ii),x)
	}


	expect_silent(fs<-makeFilterStats(cc,filterStats="var"))
	expect_silent(fs<-makeFilterStats(cc,filterStats=c("mean","var")))
	expect_silent(fs<-makeFilterStats(cc))
	expect_silent(fs<-makeFilterStats(cc,whichClusterIgnore="Cluster1"))
	
	expect_silent(fs<-makeFilterStats(se))
	expect_silent(fs<-makeFilterStats(se,filterStats="var"))
	expect_silent(fs<-makeFilterStats(se,filterStats=c("mean","var")))
	expect_equal(sort(filterNames(fs)),sort(c("b","mean","var")))
	expect_silent(fs<-makeFilterStats(fs,filterStats=c("var")))
	expect_equal(sort(filterNames(fs)),sort(c("b","mean","var")))
	expect_silent(fs<-makeFilterStats(fs,filterStats=c("abscv")))
	expect_equal(sort(filterNames(fs)),sort(c("b","mean","var","abscv")))
	
	expect_silent(fs<-makeFilterStats(sceSimData,filterStats="var"))
	expect_silent(fs<-makeFilterStats(sceSimData,filterStats=c("mean","var")))
	expect_silent(fs<-makeFilterStats(sceSimData))
	
	###Check getting out correct ones.
	contData<-simData[,1:20]
    expect_silent(temp<-makeFilterStats(contData))
	v<-apply(contData,1,var)
	expect_equal(clusterExperiment:::filterStats(temp,type="var"),DataFrame("var"=v))
	#should keep top 10
	expect_silent(x<-filterData(temp,percentile=10,filterStats="var",keepLarge=TRUE,absolute=FALSE))
	whTop<-order(v,decreasing=TRUE)[1:10]
	expect_equal(DataFrame("var"=v[whTop]),clusterExperiment:::filterStats(x,type="var"))
	expect_equal(contData[whTop,],unname(assay(x)))

	###Check getting out correct ones with transformation (isCount=TRUE)
	countData<-simCount[,1:20]
    expect_silent(temp<-makeFilterStats(countData,isCount=TRUE))
	v<-apply(log2(countData+1),1,var)
	expect_equal(clusterExperiment:::filterStats(temp,type="var"),DataFrame("var"=v))
	#should keep top 10
	expect_silent(x<-filterData(temp,percentile=10,filterStats="var",keepLarge=TRUE,absolute=FALSE))
	whTop<-order(v,decreasing=TRUE)[1:10]
	expect_equal(DataFrame("var"=v[whTop]),clusterExperiment:::filterStats(x,type="var"))
	expect_equal(log2(countData[whTop,]+1),unname(transformData(x,isCount=TRUE)))
	
	
})

test_that("defaultNDims works",{
    nDim<-3
    expect_silent(dr3<-makeReducedDims(ceSimData,reducedDims="PCA",maxDims=nDim))
    expect_equal(defaultNDims(dr3,"PCA"),3)
	
	expect_equal(defaultNDims(SingleCellExperiment(mat),"PCA"),min(dim(mat)))
	expect_equal(defaultNDims(SingleCellExperiment(ceSimData),"PCA"),50)
	
})