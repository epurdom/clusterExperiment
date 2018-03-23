context("clusterSingle")



test_that("`clusterSingle` works with matrix, ClusterExperiment objects,
          SummarizedExperiments", {
  #---
  #Matrix
  #---
	expect_silent(clustNothing <- clusterSingle(mat,
	                           subsample=FALSE, sequential=FALSE,
	                           mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"), isCount=FALSE))
	expect_equal(clusterLabels(clustNothing),"clusterSingle")
	expect_is(clustNothing, "ClusterExperiment")
	expect_is(clustNothing, "SummarizedExperiment")
	#check get the same result from previous times
	expect_equal(primaryCluster(clustNothing),c(1,2,3,1,3,1,1,1,1,2,2,1,1,1,1))

	#test clusterLabel
	expect_silent(clustNothing2 <- clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(k=3), clusterFunction="pam"), subsample=FALSE,
		sequential=FALSE, isCount=FALSE, clusterLabel="myownClustering"))
	expect_equal(clusterLabels(clustNothing2),"myownClustering")


	#test default 01 distance
	expect_silent(clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(alpha=0.1),clusterFunction="tight"),
	                              subsample=FALSE, sequential=FALSE,
	                              isCount=FALSE))
	#error because not 01 distance
	expect_error(clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(alpha=0.1),clusterFunction="tight",distFunction=function(x){dist(x,method="manhattan")}),
	                           subsample=FALSE, sequential=FALSE,isCount=FALSE),"distance function must give values between 0 and 1")

	#test default K distance
	expect_silent(clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(k=3),clusterFunction="hierarchicalK"),subsample=FALSE, sequential=FALSE, isCount=FALSE))

	#warn wrong arguments
	expect_warning(clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(k=3,alpha=0.1),clusterFunction="tight"),
	                  subsample=FALSE, sequential=FALSE,
	                  ,isCount=FALSE),"arguments passed via clusterArgs to the clustering function tight are not all applicable")
	#turn off warning
	expect_silent(clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(k=3,alpha=0.1),checkArgs=FALSE,clusterFunction="tight"),
	                  subsample=FALSE, sequential=FALSE,
	                  ,isCount=FALSE))

  #---
  #SE/SCE
  #---
	###Apply to SE
	expect_silent(clustNothing2 <- clusterSingle(se,
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
	    subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(clusterMatrix(clustNothing2), clusterMatrix(clustNothing))

	###Apply to SCE
	expect_silent(clustNothing3 <- clusterSingle(sce, mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
	    subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(clusterMatrix(clustNothing2), clusterMatrix(clustNothing3))
	expect_silent(clustNothing3 <- clusterSingle(sce,
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
		subsample=FALSE, sequential=FALSE, isCount=FALSE))

	#---
	#CE
	#---
	#test running on ClusterExperiment Object -- should add the new clustering
	expect_error(clustNothing3 <- clusterSingle(clustNothing2, mainClusterArgs=list(clusterArgs=list(k=4),clusterFunction="pam"),
	                              subsample=FALSE, sequential=FALSE,
	                              isCount=FALSE),"setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed")
	expect_silent(clustNothing3 <- clusterSingle(clustNothing2, mainClusterArgs=list(clusterArgs=list(k=4),clusterFunction="pam"),
	                              subsample=FALSE, sequential=FALSE ))
	expect_equal(NCOL(clusterMatrix(clustNothing3)),2)
	expect_equal(length(table(primaryCluster(clustNothing3))),4,info="Check reset primary cluster after run clusterSingle")

})

test_that("`clusterSingle` works with reduceMethod", {
	#existing dim reduce values
  	expect_silent(clustNothing4 <- clusterSingle(sceSimDataDimRed,
		reduceMethod="PCA", nDims=3,
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
		subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(sort(reducedDimNames(clustNothing4)),sort(c("PCA","TSNE")))
	expect_silent(clustNothing5 <- clusterSingle(t(reducedDim(sceSimDataDimRed,"PCA")[,1:3]),
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
        subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(clusterMatrix(clustNothing4), clusterMatrix(clustNothing5))

	#test nDims=NA
  	expect_silent(clustNothing1 <- clusterSingle(sceSimDataDimRed,
		reduceMethod="PCA", nDims=NA,
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
		subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(sort(reducedDimNames(clustNothing4)),sort(c("PCA","TSNE")))
	expect_silent(clustNothing2 <- clusterSingle(t(reducedDim(sceSimDataDimRed,"PCA")),
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
        subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(clusterMatrix(clustNothing1), clusterMatrix(clustNothing2))

	#blank sceSimData object
  	expect_silent(clustNothing6 <- clusterSingle(sceSimData,
		reduceMethod="PCA", nDims=3,
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
		subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(sort(reducedDimNames(clustNothing6)),sort(c("PCA")))
	expect_silent(clustNothing7 <- clusterSingle(t(reducedDim(sceSimDataDimRed,"PCA")[,1:3]),
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
        subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(clusterMatrix(clustNothing6), clusterMatrix(clustNothing7))

  	expect_silent(clustNothing8 <- clusterSingle(seSimData,
		reduceMethod="PCA", nDims=3,
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
		subsample=FALSE, sequential=FALSE, isCount=FALSE))
	expect_equal(sort(reducedDimNames(clustNothing8)),sort(c("PCA")))
	expect_equal(clusterMatrix(clustNothing8), clusterMatrix(clustNothing7))
	#check consistency
	expect_equal(primaryCluster(clustNothing8),rep(c(1,2,3),each=100))

  	expect_warning(clusterSingle(sceSimData,
		reduceMethod="PCA", nDims=NA,
		mainClusterArgs=list(clusterArgs=list(k=3),clusterFunction="pam"),
		subsample=FALSE, sequential=FALSE, isCount=FALSE),"all singular values are requested")



})

test_that("`clusterSingle` works with hdf5Matrix",{
	# clusterSingle(x, diss, subsample = TRUE,
	#   sequential = FALSE, mainClusterArgs = NULL, subsampleArgs = NULL,
	#   seqArgs = NULL, isCount = FALSE, transFun = NULL,
	#   reduceMethod = c("none", listBuiltInReducedDims(),
	#   listBuiltInFilterStats()), nDims = defaultNDims(x, reduceMethod),
	#   clusterLabel = "clusterSingle", checkDiss = TRUE)
	#
    for(kk in 1:length(listBuiltInFunctions)){
      expect_silent(clustNothing1 <- clusterSingle(sceSimDataDimRed,
  		  reduceMethod = "none",
    		  mainClusterArgs=list( clusterArgs=list(k=3),clusterFunction=listBuiltInFunctions()[[kk]]),
    	       subsample=FALSE, sequential=FALSE, isCount=FALSE))  	
		
  	  expect_silent(clustNothing2 <- clusterSingle(hdfSCE,
		  reduceMethod = "none",
  		  mainClusterArgs=list( clusterArgs=list(k=3),clusterFunction=listBuiltInFunctions()[[kk]]),
  	       subsample=FALSE, sequential=FALSE, isCount=FALSE))
		expect_equal(clusterMatrix(clustNothing1) ,clusterMatrix(clustNothing2))
    }
	
})
test_that("`clusterSingle` works with filtering", {
    ####Check built in functions ####
	for(fs in listBuiltInFilterStats()){
	    expect_silent(clusterSingle(mat,
	        subsample=FALSE, sequential=FALSE, reduceMethod=fs,
	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	 	   isCount=FALSE))

	}
    expect_warning(clusterSingle(mat,  subsample=FALSE, sequential=FALSE,
           reduceMethod="var", nDims=NROW(mat)+1,
		   mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
		 isCount=FALSE),
                 "the number of most features requested after filtering is either missing or larger than the number of features. Will not do any filtering")
    expect_warning(clusterSingle(mat,
      subsample=FALSE, sequential=FALSE,
      reduceMethod="none",nDims=3,
	  mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),isCount=FALSE),
                   "specifying nDims has no effect if reduceMethod==`none`")

    #check returns filterStats
    expect_silent(cc <- clusterSingle(mat,
    subsample=FALSE, sequential=FALSE,
    reduceMethod="var", nDims=3,
    mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), isCount=FALSE))
  	expect_equal(sort(filterNames(cc)),"var")

    expect_silent(cc2 <- clusterSingle(se,
    subsample=FALSE, sequential=FALSE,
    reduceMethod="var", nDims=3,
    mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), isCount=FALSE))
  	expect_equal(sort(filterNames(cc2)),c("b","var"))

  	expect_silent(cc3 <- clusterSingle(sceSimData,
    subsample=FALSE, sequential=FALSE,
    reduceMethod="var", nDims=3,
    mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), isCount=FALSE))
  	expect_equal(sort(filterNames(cc3)),c("b","var"))

  	#test using existing values
  	expect_silent(cc4 <- clusterSingle(sceFull,
    subsample=FALSE, sequential=FALSE,
    reduceMethod="var", nDims=3,
    mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), isCount=FALSE))
  expect_equal(sort(filterNames(cc4)),sort(c("var","b","Filter1","Filter2")))


})


test_that("`clusterSingle` consistent results (no transformation)", {

	#----
	#check get the same result from previous times when transform the data
    #----
	#make it big enough can do pca and filter...
	contData<-simData[,1:20]
    testSE<-SummarizedExperiment(contData)
    testSCE<-as(testSE,"SingleCellExperiment")
	expectTrans1<-round(contData[1,],2)

	###
	#var filter
	###

	#matrix
	expect_silent(cc<-clusterSingle(contData,
	        subsample=FALSE, sequential=FALSE, reduceMethod="var",
	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	 	   isCount=FALSE))
    expect_equal(primaryCluster(cc),c(1,2,2,1,3,3,1,3,1,2,3,3,2,1,1,3,2,2,3,2))
	expect_equal(round(transformData(cc)[1,],2), expectTrans1)

	expect_silent(cc<-clusterSingle(contData,
	        subsample=FALSE, sequential=FALSE, reduceMethod="var",
	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	 	   isCount=FALSE))
    expect_equal(primaryCluster(cc),c(1,2,2,1,3,3,1,3,1,2,3,3,2,1,1,3,2,2,3,2))


   #SE
   expect_silent(cc2<-clusterSingle(testSE,
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=FALSE))
   expect_equal(primaryCluster(cc2),primaryCluster(cc))


   #SCE
   expect_silent(cc3<-clusterSingle(testSCE,
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=FALSE))
   expect_equal(primaryCluster(cc3),primaryCluster(cc))


   #CE
   expect_silent(cc5<-clusterSingle(cc,
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),clusterLabel="redoCS"))
   expect_equal(primaryCluster(cc5),primaryCluster(cc))




	###
	#PCA
	###

  #matrix
  expect_silent(cc<-clusterSingle(contData,
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=FALSE))
  expect_equal(primaryCluster(cc),c(1,2,1,2,1,3,2,2,2,1,3,3,3,1,2,3,3,2,3,1))
	expect_equal(round(transformData(cc)[1,],2), expectTrans1)


  #SE
  expect_silent(cc2<-clusterSingle(testSE,
       subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
       nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=FALSE))
  expect_equal(primaryCluster(cc2),primaryCluster(cc))


  #SCE
  expect_silent(cc3<-clusterSingle(testSCE,
       subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
       nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=FALSE))
  expect_equal(primaryCluster(cc3),primaryCluster(cc))

  #CE
  expect_silent(cc5<-clusterSingle(cc,
       subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
       nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),clusterLabel="redoCS"))
  expect_equal(primaryCluster(cc5),primaryCluster(cc))


	###
	#Nothing
	###

	#matrix
	expect_silent(cc<-clusterSingle(contData,
	      subsample=FALSE, sequential=FALSE, reduceMethod="none",
	      nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=FALSE))
	expect_equal(primaryCluster(cc),c(1,2,3,1,2,3,1,2,1,2,3,2,3,3,1,3,3,3,2,2))
	expect_equal(round(transformData(cc)[1,],2), expectTrans1)


	#SE
	expect_silent(cc2<-clusterSingle(testSE,
	     subsample=FALSE, sequential=FALSE, reduceMethod="none",
	     nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=FALSE))
	expect_equal(primaryCluster(cc2),primaryCluster(cc))


	#SCE
	expect_silent(cc3<-clusterSingle(testSCE,
	     subsample=FALSE, sequential=FALSE, reduceMethod="none",
	     nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=FALSE))
	expect_equal(primaryCluster(cc3),primaryCluster(cc))


	#CE
	expect_silent(cc5<-clusterSingle(cc,
	     subsample=FALSE, sequential=FALSE, reduceMethod="none",
	     nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),clusterLabel="redoCS"))
	expect_equal(primaryCluster(cc5),primaryCluster(cc))

})


test_that("`clusterSingle` consistent results (with transformation)", {

	#----
	#check get the same result from previous times when transform the data
    #----
	#make it big enough can do pca and filter...
	countData<-simCount[,1:20]
    testSE<-SummarizedExperiment(countData)
    testSCE<-as(testSE,"SingleCellExperiment")
	expectTrans1<-round(log2(countData[1,]+1),2)
	###
	#var filter
	###

	#matrix
	expect_silent(cc<-clusterSingle(countData,
	        subsample=FALSE, sequential=FALSE, reduceMethod="var",
	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	 	   isCount=TRUE))
   expect_equal(primaryCluster(cc),c(1,1,2,1,2,1,2,2,3,2,1,3,1,3,1,1,2,3,2,2))
	expect_equal(round(transformData(cc)[1,],2), expectTrans1)

   #SE
   expect_silent(cc2<-clusterSingle(testSE,
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
   expect_equal(primaryCluster(cc2),primaryCluster(cc))


   #SCE
   expect_silent(cc3<-clusterSingle(testSCE,
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
   expect_equal(primaryCluster(cc3),primaryCluster(cc))

   #CE
   expect_silent(cc5<-clusterSingle(cc,
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),clusterLabel="redoCS"))
   expect_equal(primaryCluster(cc5),primaryCluster(cc))




	###
	#PCA
	###

  #matrix
  expect_silent(cc<-clusterSingle(countData,
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
  expect_equal(primaryCluster(cc),c(1,2,3,3,1,3,1,3,3,3,3,2,3,2,3,3,3,2,2,3))
expect_equal(round(transformData(cc)[1,],2), expectTrans1)

  #SE
  expect_silent(cc2<-clusterSingle(testSE,
       subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
       nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=TRUE))
  expect_equal(primaryCluster(cc2),primaryCluster(cc))


  #SCE
  expect_silent(cc3<-clusterSingle(testSCE,
       subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
       nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=TRUE))
  expect_equal(primaryCluster(cc3),primaryCluster(cc))

  #CE
  expect_silent(cc5<-clusterSingle(cc,
       subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
       nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),clusterLabel="redoCS"))
  expect_equal(primaryCluster(cc5),primaryCluster(cc))


	###
	#Nothing
	###

	#matrix
	expect_silent(cc<-clusterSingle(countData,
	      subsample=FALSE, sequential=FALSE, reduceMethod="none",
	      nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=TRUE))
	expect_equal(primaryCluster(cc),c(1,2,1,2,1,2,1,2,3,1,2,3,2,3,2,2,2,3,3,1))

	#SE
	expect_silent(cc2<-clusterSingle(testSE,
	     subsample=FALSE, sequential=FALSE, reduceMethod="none",
	     nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=TRUE))
	expect_equal(primaryCluster(cc2),primaryCluster(cc))


	#SCE
	expect_silent(cc3<-clusterSingle(testSCE,
	     subsample=FALSE, sequential=FALSE, reduceMethod="none",
	     nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	   isCount=TRUE))
	expect_equal(primaryCluster(cc3),primaryCluster(cc))


	#CE
	expect_silent(cc5<-clusterSingle(cc,
	     subsample=FALSE, sequential=FALSE, reduceMethod="none",
	     nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),clusterLabel="redoCS"))
	expect_equal(primaryCluster(cc5),primaryCluster(cc))

})
		  # > clustSeqHier_v2 <- clusterSingle(simData,
		  # + sequential=FALSE, subsample=TRUE, subsampleArgs=list(resamp.n=100, samp.p=0.7,
		  # + clusterFunction="kmeans", clusterArgs=list(nstart=10)),
		  # + seqArgs=list(beta=0.8, k0=5), mainClusterArgs=list(minSize=5,clusterFunction="hierarchical01"))
		  # Error in .local(x, diss, ...) :
		  #   For the clusterFunction algorithm type (' 01 ') given in 'mainClusterArgs', must supply arguments: alpha These must be supplied as elements of the list of 'clusterArgs' given in 'mainClusterArgs'
		  # > set.seed(44261)
		  # > clustSeqHier_v2 <- clusterSingle(simData,
		  # + sequential=FALSE, subsample=TRUE, subsampleArgs=list(resamp.n=100, samp.p=0.7,
		  # + clusterFunction="kmeans", clusterArgs=list(nstart=10)),
		  # + seqArgs=list(beta=0.8, k0=5), mainClusterArgs=list(minSize=5,clusterFunction="hierarchical01",clusterArgs=list(alpha=0.1)))


test_that("Different options algorithms of `mainClustering` ", {
  #check builtIn algorithms
  #bigger matrix so not kill spectral
  set.seed(3325)
  biggerMat<-matrix(data=rnorm(20*50), ncol=50)

  kMethods<-listBuiltInTypeK()
	for(cf in kMethods){
	    expect_silent(clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(k=3), clusterFunction=cf),
	  			subsample=FALSE, sequential=FALSE,isCount=FALSE)
				)
		#post-processing arguments for type 'K'
		#Upped
		expect_silent(clusterSingle(biggerMat, mainClusterArgs= list(clusterArgs=list(k=3), clusterFunction=cf,findBestK=TRUE,removeSil=TRUE), subsample=FALSE, sequential=FALSE,isCount=FALSE))

	  }
      aMethods<-listBuiltInType01()
  	for(cf in aMethods){
  	  expect_silent(clusterSingle(mat, mainClusterArgs= list(clusterArgs=list(alpha=0.1),clusterFunction=cf),
                     subsample=FALSE, sequential=FALSE,isCount=FALSE))
	 }

  ########
  #Check mainClustering
  ########
  ###Check pam exactly same:
  expect_silent(x<-mainClustering(mat, clusterFunction="pam",clusterArgs=list(k=3),
           minSize=1, removeSil=FALSE))
  expect_equal(length(x),ncol(mat))
  x2<-cluster::pam(t(mat),k=3,cluster.only=TRUE)
  expect_equal(x,x2)
  ###Check hierarchicalK exactly same:
  expect_silent(x<-mainClustering(mat, clusterFunction="hierarchicalK",clusterArgs=list(k=3),
              minSize=1, removeSil=FALSE))
  expect_equal(length(x),ncol(mat))
  x2<-stats::cutree(stats::hclust(dist(t(mat))),k=3)
  expect_equal(x,x2)


  #check giving wrong parameters gives warning:
  expect_warning(mainClustering(mat, clusterFunction="tight", clusterArgs=list(alpha=0.1),
      minSize=5, removeSil=TRUE),"do not match the algorithmType")
  expect_warning(mainClustering(mat, clusterFunction="tight", clusterArgs=list(alpha=0.1,k=3),
      minSize=5, removeSil=TRUE),"do not match the algorithmType")
  expect_error(mainClustering(mat, clusterFunction="tight",
	      minSize=5),"must supply arguments alpha")
  expect_error(mainClustering(mat, clusterFunction="pam", clusterArgs=list(alpha=0.1),
			      minSize=5, removeSil=TRUE),"must supply arguments k")
  expect_warning(mainClustering(mat, clusterFunction="tight", clusterArgs=list(k=3,alpha=0.1),
						      minSize=5),"arguments passed via clusterArgs to the clustering function tight are not all applicable")
  expect_warning(mainClustering(mat, clusterFunction="pam", clusterArgs=list(k=3,alpha=0.1),
				      minSize=5, removeSil=TRUE),"arguments passed via clusterArgs to the clustering function pam are not all applicable")



  expect_warning(mainClustering(mat, clusterFunction="tight", clusterArgs=list(alpha=0.1, evalClusterMethod="average")),"arguments passed via clusterArgs to the clustering function tight are not all applicable")
  expect_warning(mainClustering(mat, clusterFunction="hierarchical01", clusterArgs=list(alpha=0.1, minSize.core=4)),"arguments passed via clusterArgs to the clustering function hclust are not all applicable")

  #test default 01 distance
  expect_silent(mainClustering(mat, clusterFunction="tight", clusterArgs=list(alpha=0.1)))
  #test default K distance
  expect_silent(mainClustering(mat, clusterFunction="hierarchicalK", clusterArgs=list(k=3)))


  #check turn off if checkArgs=TRUE
  expect_silent(mainClustering(mat, clusterFunction="tight", clusterArgs=list(alpha=0.1),checkArgs=FALSE,
                          minSize=5, removeSil=TRUE))
  expect_silent(mainClustering(mat, clusterFunction="pam", clusterArgs=list(alpha=0.1),checkArgs=FALSE,
                          minSize=5, removeSil=TRUE, findBestK=TRUE))
  expect_silent(mainClustering(mat, clusterFunction="tight", clusterArgs=list(alpha=0.1,evalClusterMethod="average"),checkArgs=FALSE))
  expect_silent(mainClustering(mat, clusterFunction="hierarchical01", checkArgs=FALSE,
                          clusterArgs=list(alpha=0.1,minSize.core=4)))

})

test_that("Different options of mainClustering",{
    #check errors and warnings
    expect_error(clusterSingle(mat,  subsample=FALSE, sequential=TRUE, seqArgs=list(verbose=FALSE), isCount=FALSE,mainClusterArgs=list(clusterFunction="pam")),
                 "seqArgs must contain element 'k0'")
    expect_error(clusterSingle(mat,  subsample=FALSE, sequential=TRUE, seqArgs=list(verbose=FALSE), isCount=FALSE, mainClusterArgs=list(clusterFunction="pam","findBestK"==TRUE)),
                 "seqArgs must contain element 'k0'")
    expect_error(clusterSingle(mat,
                                 subsample=FALSE, sequential=FALSE,
                                 mainClusterArgs=list(clusterFunction="tight",clusterArgs=list(k=3)), isCount=FALSE),
                   "must supply arguments: alpha")
    expect_warning(clusterSingle(mat,  subsample=FALSE, sequential=FALSE, mainClusterArgs=list(clusterFunction="tight",clusterArgs=list(alpha=0.1),findBestK=TRUE),isCount=FALSE),
                   "Some arguments passed via '...' in mainClustering do not match the algorithmType")
    expect_error(clusterSingle(mat,
                                 subsample=FALSE, sequential=FALSE,
                                 mainClusterArgs=list(clusterFunction="tight",clusterArgs=list(alpha=0.1),distFunction=function(x){abs(cor(t(x)))}),isCount=FALSE),
                   "Dissimilarity matrix must have zero values on the diagonal")


})

test_that("Different options of subsampling",{
	expect_warning(clustSubsample <- clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=3, clusterArgs=list(k=3)), mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)), isCount=FALSE),"a clusterFunction was not set for subsampleClustering")
    expect_equal(NCOL(coClustering(clustSubsample)),NCOL(mat))

    #check subsample works with all of the builtin functions and opposite type in mainClusterArgs
    set.seed(3325)
    biggerMat<-matrix(data=rnorm(20*100), ncol=100)
    kMethods<-listBuiltInTypeK()
	for(cf in kMethods){
		set.seed(1045)
	    expect_silent(clusterSingle(biggerMat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=20, clusterArgs=list(k=3),clusterFunction=cf,classifyMethod="InSample"), mainClusterArgs=list(clusterFunction="hierarchical01", clusterArgs=list(alpha=0.3)),isCount=FALSE))
       if(!is.null(getBuiltInFunction(cf)@classifyFUN)){
   		set.seed(1045)
	    expect_silent(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=20, clusterArgs=list(k=3),clusterFunction=cf,classifyMethod="All"), mainClusterArgs=list(clusterFunction="hierarchical01", clusterArgs=list(alpha=0.1)),isCount=FALSE))
		set.seed(1045)
	    expect_silent(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=40, clusterArgs=list(k=3),clusterFunction=cf,classifyMethod="OutOfSample"), mainClusterArgs=list(clusterFunction="hierarchical01", clusterArgs=list(alpha=0.1)),isCount=FALSE))

       }
	}
    aMethods<-listBuiltInType01()
	for(cf in aMethods){

		set.seed(1045)
	    expect_silent(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=20, clusterArgs=list(alpha=0.1),clusterFunction=cf,classifyMethod="InSample"), mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)),isCount=FALSE))
        if(!is.null(getBuiltInFunction(cf)@classifyFUN)){
			##Check outofsample/all
			set.seed(1045)
	    	expect_silent(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=20, clusterArgs=list(alpha=0.1),clusterFunction=cf,classifyMethod="All"), mainClusterArgs=list(clusterFunction="hierarchical01", clusterArgs=list(k=3)),isCount=FALSE))
			set.seed(1045)
	    	expect_silent(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=40, clusterArgs=list(alpha=0.1),clusterFunction=cf,classifyMethod="OutOfSample"), mainClusterArgs=list(clusterFunction="hierarchical01", clusterArgs=list(k=3)),isCount=FALSE))

        }
	}

    ## get NA values
	set.seed(1045)
	# with the new version, we fix NA's in subsampleClustering
#     expect_error(clusterSingle(mat,
#        subsample=TRUE, sequential=FALSE,
#        subsampleArgs=list(resamp.num=20,clusterArgs=list(k=3),clusterFunction="pam",
# 	   classifyMethod="OutOfSample"),
#        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),isCount=FALSE),"NA values found in dissimilarity matrix")

    #warnings in missing args in subsample -- should borrow from mainClusterArgs .
    expect_warning(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(clusterFunction="pam",resamp.num=3),  mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), isCount=FALSE),
                   "missing arguments k provided from those in 'mainClusterArgs'")
	#warnings in missing clusterFunction in subsample -- should borrow from mainClusterArgs .
	expect_warning(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(resamp.num=3,clusterArgs=list(k=3)),  mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), isCount=FALSE),
			                      "a clusterFunction was not set for subsampleClustering")
	#different function types -- should error out.
    expect_error(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(clusterFunction="pam",resamp.num=3),  mainClusterArgs=list(clusterFunction="tight",clusterArgs=list(alpha=0.1)), isCount=FALSE),
	"must supply arguments: k")


})


test_that("Different passed options of seqCluster",{
    #check sequential
    expect_silent(clustSeq <- clusterSingle(mat,subsample=FALSE, sequential=TRUE,mainClusterArgs=list(clusterFunction="pam"),isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)))
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,mainClusterArgs=list(clusterFunction="pam"),isCount=FALSE), "if sequential=TRUE, must give seqArgs so as to identify k0 and beta")
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,mainClusterArgs=list(clusterFunction="pam"),isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE)), "seqArgs must contain element 'beta'")
	expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,mainClusterArgs=list(clusterFunction="pam"),isCount=FALSE,seqArgs=list(beta=0.9,verbose=FALSE)), "seqArgs must contain element 'k0'")

	#right clusterFunctions
	expect_error(clusterSingle(mat, mainClusterArgs=list(clusterFunction="kmeans"), subsample=TRUE, sequential=TRUE, subsampleArgs=list(clusterFunction="pam",n.sample=40), isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
							  "If choosing subsample=TRUE, the clusterFunction used in the mainClustering step must take input that is dissimilarity")
	expect_error(clusterSingle(mat, mainClusterArgs=list(clusterFunction="tight"),  subsample=FALSE, sequential=TRUE, isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
		"if subsample=FALSE, sequentical clustering can only be implemented with a clusterFunction with algorithmType 'K'")
	#warning if try to set k
	expect_warning(clusterSingle(mat, mainClusterArgs=list(clusterFunction="pam"), subsample=TRUE, sequential=TRUE, subsampleArgs=list(clusterFunction="pam",n.sample=40,clusterArgs=list(k=3)), isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
								  "Setting 'k' in subsampleArgs when sequential=TRUE is called will have no effect.")
	expect_warning(clusterSingle(mat, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), subsample=FALSE, sequential=TRUE, subsampleArgs=list(clusterFunction="pam",n.sample=40), isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
							  								  "Setting 'k' in mainClusterArgs when sequential clustering is requested will have no effect.")

	#check all algorithms
    kMethods<-listBuiltInTypeK()
	for(cf in kMethods){
		#check if no subsampling
		expect_silent(clusterSingle(mat, mainClusterArgs=list(clusterFunction=cf),
	                              subsample=FALSE, sequential=TRUE,
	                              isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)))
	}
	kMethods<-listBuiltInType01()
	for(cf in kMethods){
		#check if no subsampling
		expect_silent(clusterSingle(mat, mainClusterArgs=list(clusterFunction=cf),
	                              subsample=TRUE, sequential=TRUE,
								  subsampleArgs=list(clusterFunction="pam",n.sample=40),
	                              isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)))
	}


})

test_that("Different direct options of `clusterSingle` ", {
  #check isCount
  expect_silent(clusterSingle(smSimCount,
                           subsample=FALSE, sequential=FALSE,
                           mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"),isCount=TRUE) )
  expect_error(clusterSingle(smSimData,
                          subsample=FALSE, sequential=FALSE,
                          mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"),isCount=TRUE),"User-supplied `transFun` produces NA values",info="test error handling for isCount=TRUE when can't take log")


})

test_that("`clusterSingle` preserves the colData and rowData of SE", {
  expect_silent(cl<-clusterSingle(se,
                      subsample=FALSE, sequential=FALSE,
                      mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),isCount=FALSE))

  expect_equal(colData(cl),colData(se))
  expect_equal(rownames(cl),rownames(se))
  expect_equal(colnames(cl),colnames(se))
  expect_equal(metadata(cl),metadata(se))
  expect_equal(rowData(cl),rowData(se))

})

## Windows does not support mclapply
if(.Platform$OS.type == "unix"){
  test_that("`clusterSingle` works with parallel subsampling", {
    expect_silent(clustSubsample2 <- clusterSingle(mat,  subsample=TRUE,
                                                   sequential=FALSE,
                                                   subsampleArgs=list(clusterFunction="kmeans",
                                                                      resamp.num=3,
                                                                      clusterArgs=list(k=3),
                                                                      ncores=2),
                                                   mainClusterArgs=list(clusterFunction="pam",
                                                                        clusterArgs=list(k=3)),
                                                   isCount=FALSE))

    expect_silent(clustSubsample1 <- clusterSingle(mat,  subsample=TRUE,
                                                   sequential=FALSE,
                                                   subsampleArgs=list(clusterFunction="kmeans",
                                                                      resamp.num=3,
                                                                      clusterArgs=list(k=3),
                                                                      ncores=1),
                                                   mainClusterArgs=list(clusterFunction="pam",
                                                                        clusterArgs=list(k=3)),
                                                   isCount=FALSE))
  })
}
