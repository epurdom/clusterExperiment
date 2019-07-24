context("clusterSingle")



   
test_that("`clusterSingle` works with matrix, ClusterExperiment objects,
          SummarizedExperiments", {
    #---
    #Matrix
    #---
	expect_silent(clustNothing <- 
        clusterSingle(mat, inputType="X",
        subsample=FALSE, sequential=FALSE,
        mainClusterArgs=list(clusterArgs=list(k=3), 
            clusterFunction="pam"), 
        isCount=FALSE))
	expect_equal(clusterLabels(clustNothing),"clusterSingle")
	expect_is(clustNothing, "ClusterExperiment")
	expect_is(clustNothing, "SummarizedExperiment")
	truth<- unname(cluster::pam(x=t(mat),cluster.only=TRUE,k=3))
    #in case change in sim, make sure same as from previous tests
    expect_equal(truth,c(1,2,3,1,3,1,1,1,1,2,2,1,1,1,1)) 
	expect_equal(primaryCluster(clustNothing),truth)

	#test clusterLabel
	expect_silent(clustNothing2 <- clusterSingle(mat,
        inputType="X",
        mainClusterArgs= list(clusterArgs=list(k=3), clusterFunction="pam"), 
        subsample=FALSE,
		sequential=FALSE, isCount=FALSE, 
        clusterLabel="myownClustering"))
	expect_equal(clusterLabels(clustNothing2),"myownClustering")


	#test default 01 distance
	expect_silent(clusterSingle(inputMatrix=dissMat, inputType="diss", mainClusterArgs= list(clusterArgs=list(alpha=0.1),clusterFunction="tight"),
	                              subsample=FALSE, sequential=FALSE,
	                              isCount=FALSE))

	#test default K distance
	expect_silent(clusterSingle(dissMat, inputType="diss", mainClusterArgs= list(clusterArgs=list(k=3),clusterFunction="hierarchicalK"),subsample=FALSE, sequential=FALSE, isCount=FALSE))

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

test_that("`clusterSingle` warnings/errors work", {    
	expect_error(clusterSingle(mat, inputType="X", 
        mainClusterArgs= list(clusterArgs=list(alpha=0.1),
            clusterFunction="tight"),
        distFunction=function(x){dist(x,method="manhattan")},
        subsample=FALSE, checkDiss=TRUE,
        sequential=FALSE,isCount=FALSE),
    "distance function must give values between 0 and 1")
	#test warn with wrong arguments
	expect_warning(clusterSingle(dissMat, inputType="diss", mainClusterArgs= list(clusterArgs=list(k=3,alpha=0.1),clusterFunction="tight"),
	                  subsample=FALSE, sequential=FALSE,
	                  ,isCount=FALSE),"arguments passed via clusterArgs to the clustering function tight are not all applicable")
	#test turning off warning
	expect_silent(clusterSingle(dissMat, 
        inputType="diss", 
        mainClusterArgs= list(clusterArgs=list(k=3,alpha=0.1), clusterFunction="tight"),
        warnings=FALSE,
	    subsample=FALSE, 
        sequential=FALSE,
        isCount=FALSE))    
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
    ## Loops over all built-in clusterFunctions and makes sure they work with hdf5
    kMethods<-listBuiltInTypeK()
    aMethods<-listBuiltInType01()
	seedValue<-571839
    kClusterArgs<-list(k=3)
    aClusterArgs<-list(alpha=0.3)
  	for(cf in c(aMethods,kMethods)){
		#print(cf)
        clArgs<-switch(algorithmType(cf),
            "K"=kClusterArgs,
            "01"=aClusterArgs)
        ### Annoyingly, have to separately for "X" and "diss"
		if("X" %in% inputType(cf)){
    		set.seed(seedValue)
  	        expect_silent(clust1<-clusterSingle(sceSimDataDimRed,   
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs,
                    clusterFunction=cf),
      	  		subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE)
      		)
    		set.seed(seedValue)
      	    expect_silent(clust2<-clusterSingle(hdfSCE, 
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs,
                    clusterFunction=cf),
                subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE)
      		)
    		set.seed(seedValue)
      	    expect_silent(clust3<-clusterSingle(hdfObj, 
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs, 
                    clusterFunction=cf),
                subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE)
      		)
    		set.seed(seedValue)
      	    expect_silent(clust4<-clusterSingle(assay(hdfSCE), 
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs, 
                    clusterFunction=cf),
                subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE)
      		)
  	    }
        else{
    		set.seed(seedValue)
  	        expect_message(clust1<-clusterSingle(sceSimDataDimRed,   
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs, 
                    clusterFunction=cf),
      	  		subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE),"Note: Making nxn dissimilarity matrix."
      		)
    		set.seed(seedValue)
  	        expect_message(clust2<-clusterSingle(hdfSCE, 
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs, 
                    clusterFunction=cf),
                subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE),"Note: Making nxn dissimilarity matrix."
      		)
			
    		set.seed(seedValue)
      	    expect_message(clust3<-clusterSingle(hdfObj, 
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs, 
                    clusterFunction=cf),
                subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE),"Note: Making nxn dissimilarity matrix."
      		)
    		set.seed(seedValue)
      	    expect_message(clust4<-clusterSingle(assay(hdfSCE), 
                reduceMethod = "none", 
                mainClusterArgs= list(clusterArgs=clArgs, 
                    clusterFunction=cf),
                subsample=FALSE, 
                makeMissingDiss=TRUE,
                sequential=FALSE,
                isCount=FALSE),"Note: Making nxn dissimilarity matrix."
      		)
        }
		expect_equal(clusterMatrix(clust1) ,clusterMatrix(clust2))
		expect_equal(clusterMatrix(clust1) ,clusterMatrix(clust3))
		expect_equal(clusterMatrix(clust1) ,clusterMatrix(clust4))
  	}
 	 ####Test sequential option,
     expect_silent(clustSeq <- clusterSingle(hdfObj,
         reduceMethod="none",
         subsample=FALSE, 
         sequential=TRUE,
         mainClusterArgs=list(clusterFunction="pam"),
         isCount=FALSE,
         seqArgs=list(k0=5,beta=0.9,verbose=FALSE)
     ))
     expect_silent(clustSeq <- clusterSingle(assay(hdfObj),
         reduceMethod="none",
         subsample=FALSE, sequential=TRUE,
         mainClusterArgs=list(clusterFunction="pam"),
         isCount=FALSE,
         seqArgs=list(k0=5,beta=0.9,verbose=FALSE)))


	 ####Test subsample option
     expect_silent(clusterSingle(hdfObj, 
        reduceMethod="none", 
	 	subsample=TRUE, 
        sequential=FALSE,
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
		isCount=FALSE,
	 	subsampleArgs=list(clusterFunction="pam",resamp.num=3, 
            clusterArgs=list(k=3))
		)
	 )
	 expect_silent(clusterSingle(assay(hdfObj), reduceMethod="none", 
	 	subsample=TRUE, sequential=FALSE,
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
		isCount=FALSE,
	 	subsampleArgs=list(clusterFunction="pam",resamp.num=3, clusterArgs=list(k=3))
		)
	 )
	
})

test_that("`clusterSingle` works with filtering", {
    ####Check built in functions ####
	for(fs in listBuiltInFilterStats()){
	    expect_silent(clusterSingle(mat,
	        subsample=FALSE, 
            sequential=FALSE, 
            reduceMethod=fs,
	        nDims=3,
             mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)),
            isCount=FALSE))

	}
    expect_warning(clusterSingle(mat,  
        subsample=FALSE, 
        sequential=FALSE,
        reduceMethod="var",
        nDims=NROW(mat)+1,
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
		isCount=FALSE),
        "the number of most features requested after filtering is either missing or larger than the number of features. Will not do any filtering")
    expect_warning(clusterSingle(mat,
      subsample=FALSE, sequential=FALSE,
      reduceMethod="none",nDims=3,
	  mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
      isCount=FALSE),
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
    expect_equal(primaryCluster(cc),c(3,1,1,3,2,2,3,2,3,1,2,2,1,3,3,2,1,1,2,1))
	expect_equal(round(transformData(cc)[1,],2), expectTrans1)

	expect_silent(cc<-clusterSingle(contData,
	        subsample=FALSE, sequential=FALSE, reduceMethod="var",
	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	 	   isCount=FALSE))
    expect_equal(primaryCluster(cc),c(3,1,1,3,2,2,3,2,3,1,2,2,1,3,3,2,1,1,2,1))


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
  expect_equal(primaryCluster(cc),c(3,1,3,1,3,2,1,1,1,3,2,2,2,3,1,2,2,1,2,3))
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
	expect_equal(primaryCluster(cc),c(3,2,1,3,2,1,3,2,3,2,1,2,1,1,3,1,1,1,2,2))
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
  expect_equal(primaryCluster(cc),c(3,2,1,1,3,1,3,1,1,1,1,2,1,2,1,1,1,2,2,1))
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
	expect_equal(primaryCluster(cc),c(2,1,2,1,2,1,2,1,3,2,1,3,1,3,1,1,1,3,3,2))

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


test_that("Different clustering algorithms of `mainClustering` ", {
    #check builtIn algorithms
    #bigger matrix so not kill spectral
    #FIXME: Is this duplicative of previous test in hdf5?
    set.seed(3325)
    biggerMat<-matrix(data=rnorm(20*50), ncol=50)

    kMethods<-listBuiltInTypeK()
    for(cf in kMethods){
        if("X" %in% inputType(cf)){
            expect_silent(clusterSingle(mat,
                  makeMissingDiss=TRUE, 
                  mainClusterArgs= list(clusterArgs=list(k=3), clusterFunction=cf),
                  subsample=FALSE, sequential=FALSE,isCount=FALSE)
            )
            #add post-processing arguments for type 'K'
            #Upped
            expect_silent(clusterSingle(biggerMat, 
                  makeMissingDiss=TRUE,
                  mainClusterArgs= list(clusterArgs=list(k=3),
                      clusterFunction=cf,findBestK=TRUE,removeSil=TRUE), 
                  subsample=FALSE, sequential=FALSE,isCount=FALSE))
        }
        else{
            expect_message(clusterSingle(mat,
                  makeMissingDiss=TRUE, 
                  mainClusterArgs= list(clusterArgs=list(k=3), clusterFunction=cf),
                  subsample=FALSE, sequential=FALSE,isCount=FALSE),
                  "Making nxn dissimilarity matrix"
            )
            #add post-processing arguments for type 'K'
            #Upped
            expect_message(clusterSingle(biggerMat, 
                  makeMissingDiss=TRUE,
                  mainClusterArgs= list(clusterArgs=list(k=3), clusterFunction=cf,findBestK=TRUE,removeSil=TRUE), 
                  subsample=FALSE, sequential=FALSE,isCount=FALSE),
                  "Making nxn dissimilarity matrix"
            )
        }

    }
    aMethods<-listBuiltInType01()
    for(cf in aMethods){
        if("X" %in% inputType(cf))
            expect_silent(clusterSingle(mat, 
                makeMissingDiss=TRUE,
                mainClusterArgs= list(clusterArgs=list(alpha=0.1), clusterFunction=cf),
                subsample=FALSE, sequential=FALSE,isCount=FALSE))
        else expect_message(clusterSingle(mat, 
            makeMissingDiss=TRUE,
            mainClusterArgs= list(clusterArgs=list(alpha=0.1), clusterFunction=cf),
            subsample=FALSE, sequential=FALSE,isCount=FALSE),
            "Making nxn dissimilarity matrix")
    }
})
             
test_that("Different options of mainClustering",{
    #check errors and warnings
    expect_error(clusterSingle(mat,  
        subsample=FALSE, sequential=TRUE, 
        seqArgs=list(verbose=FALSE), 
        isCount=FALSE,
        mainClusterArgs=list(clusterFunction="pam")),
                 "required argument 'k0' is missing for the sequential clustering step")
    expect_error(clusterSingle(mat,  
        subsample=FALSE, sequential=TRUE, 
        seqArgs=list(verbose=FALSE), isCount=FALSE, 
        mainClusterArgs=list(clusterFunction="pam","findBestK"==TRUE)),
                 "required argument 'k0' is missing for the sequential clustering step")
    expect_error(clusterSingle(mat, makeMissingDiss=TRUE,
        subsample=FALSE, sequential=FALSE,
        mainClusterArgs=list(clusterFunction="tight",clusterArgs=list(k=3)), 
        isCount=FALSE),
                   "must supply arguments: alpha")
    expect_warning(clusterSingle(mat,  makeMissingDiss=TRUE,
        subsample=FALSE, sequential=FALSE, 
        mainClusterArgs=list(clusterFunction="tight",
             clusterArgs=list(alpha=0.1),findBestK=TRUE),
        isCount=FALSE),
        "Some arguments passed via mainClusterArgs in mainClustering step do not match the algorithmType of the given ClusterFunction object")
    expect_error(clusterSingle(mat,
        subsample=FALSE, sequential=FALSE,
        distFunction=function(x){abs(cor(t(x)))},  
        mainClusterArgs=list(clusterFunction="tight", clusterArgs=list(alpha=0.1)),
        checkDiss=TRUE,isCount=FALSE),
        "Dissimilarity matrix must have zero values on the diagonal")
#nReducedDims=5,k=NA,findBestK=TRUE

    #-------
    ## Check findBestK=TRUE (i.e. post-processing)
    #-------
    # Error, because needs calculate diss
    expect_error(clusterSingle(simData,
        subsample=FALSE,
        mainClusterArgs=list(clusterFunction="pam",
            kRange=2:4,
            findBestK=TRUE)),
        "Cannot do requested post processing (e.g. from arguments 'findBestK' or 'removeSil') without calculating distance matrix")
    #Should find best K is K=4
    expect_message(out<-clusterSingle(simData,
        subsample=FALSE, 
        makeMissingDiss=TRUE,
        mainClusterArgs=list(clusterFunction="pam",
            kRange=2:4,
            findBestK=TRUE)),
        "Making nxn dissimilarity matrix in order to do post-processing")
    expect_equal(length(tableClusters(out)),4)
    
    ## Should be same as if gave distance, since its pam at main clustering 
    ## (and if give distance, shouldn't get message)
    d<-as.matrix(dist(t(simData)))
    expect_silent(out2<-clusterSingle(d, inputType="diss",
            subsample=FALSE,
            mainClusterArgs=list(clusterFunction="pam",
                kRange=2:4,
                findBestK=TRUE)))
    expect_equal(primaryCluster(out),out2$clustering)
    
    ## Also same as if gave distance to mainClusterArgs
    expect_silent(out3<-clusterSingle(simData, 
            subsample=FALSE,
            mainClusterArgs=list(clusterFunction="pam",
                kRange=2:4, diss=d,
                findBestK=TRUE)))
    expect_equal(clusterMatrix(out),clusterMatrix(out3))
    
    
})

test_that("Different options of subsampling",{
    #Check saving of subsampling matrix
	expect_silent(clustSubsample <- clusterSingle(mat,  
        inputType="X",
        subsample=TRUE, 
        saveSubsamplingMatrix=TRUE,
        sequential=FALSE, makeMissingDiss=FALSE, warnings=TRUE,
        subsampleArgs=list(resamp.num=3, 
             clusterFunction="kmeans",clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="hierarchical01",
             clusterArgs=list(alpha=.3)), 
        isCount=FALSE))
    expect_equal(NROW(coClustering(clustSubsample)),nSamples(clustSubsample))
    expect_false(is.null(coClustering(clustSubsample)))
    expect_message(clustSubsample <- clusterSingle(mat,  
        subsample=TRUE, 
        saveSubsamplingMatrix=FALSE,
        sequential=FALSE, 
        subsampleArgs=list(resamp.num=3, clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)), 
        isCount=FALSE,warnings=TRUE),
        "clusterFunction was not given for subsampleClustering -- set to be the same as the mainClustering step")
    expect_true(is.null(coClustering(clustSubsample)))
    
    expect_silent(test <- clusterSingle(sce,  
        subsample=TRUE, 
        saveSubsamplingMatrix=TRUE,
        sequential=FALSE, 
        subsampleArgs=list(resamp.num=3, clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)), 
        isCount=FALSE,warnings=FALSE))
    expect_false(is.null(coClustering(test)))
    expect_silent(test <- clusterSingle(sce,  
        subsample=TRUE, 
        saveSubsamplingMatrix=FALSE,
        sequential=FALSE, 
        subsampleArgs=list(resamp.num=3, clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)), 
        isCount=FALSE,warnings=FALSE))
    expect_true(is.null(coClustering(test)))
    
    expect_silent(clustSubsampleCE <- clusterSingle(cc,  
        subsample=TRUE, 
        saveSubsamplingMatrix=TRUE,
        sequential=FALSE, 
        subsampleArgs=list(resamp.num=3, clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)),
        warnings=FALSE))
    expect_false(is.null(coClustering(clustSubsampleCE)))
    expect_silent(test <- clusterSingle(cc,  
        subsample=TRUE, 
        saveSubsamplingMatrix=FALSE,
        sequential=FALSE, 
        subsampleArgs=list(resamp.num=3, clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)),
        warnings=FALSE))
    expect_true(is.null(coClustering(test)))
    
    expect_silent(test <- clusterSingle(clustSubsampleCE,  
        subsample=TRUE, 
        sequential=FALSE, 
        saveSubsamplingMatrix=FALSE,
        subsampleArgs=list(resamp.num=3, clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)),
        warnings=FALSE))
    expect_equal(coClustering(test),coClustering(clustSubsampleCE))

	###########
    #check subsample works with all of the builtin functions and opposite type in mainClusterArgs
	###########
    set.seed(3325)
    biggerMat<-matrix(data=rnorm(20*100), ncol=100)
	#----
	# Check K methods
	#----
    kMethods<-listBuiltInTypeK()
	mainArgs<-list(clusterFunction="hierarchical01", clusterArgs=list(alpha=0.3))
	subArgs<-list(resamp.num=20, clusterArgs=list(k=3))
	for(cf in kMethods){
		set.seed(1045)
		#In sample
	    if(cf!="hierarchicalK") 
			expect_silent(clusterSingle(biggerMat,  
                subsample=TRUE, 
                sequential=FALSE, 
                makeMissingDiss=TRUE,
                subsampleArgs=c(list(clusterFunction=cf,classifyMethod="InSample"),subArgs), 
                mainClusterArgs=mainArgs,isCount=FALSE))
		else expect_message(clusterSingle(biggerMat,  
            subsample=TRUE, 
            sequential=FALSE, 
            makeMissingDiss=TRUE,
            subsampleArgs=c(list(clusterFunction=cf,classifyMethod="InSample"),subArgs), 
            mainClusterArgs=mainArgs),"Making nxn dissimilarity matrix")
       if(!is.null(getBuiltInFunction(cf)@classifyFUN)){
		#All 
   		set.seed(1045)
	    expect_silent(clusterSingle(biggerMat,  subsample=TRUE, sequential=FALSE, subsampleArgs=c(list(clusterFunction=cf,classifyMethod="All"),subArgs), mainClusterArgs=mainArgs,isCount=FALSE))
		#OutOfSample
		set.seed(1045)
	    expect_silent(clusterSingle(biggerMat,  subsample=TRUE, sequential=FALSE, subsampleArgs=c(list(clusterFunction=cf,classifyMethod="OutOfSample"),subArgs), mainClusterArgs=mainArgs,isCount=FALSE))

       }
	}
	#----
	# Check 01 methods
	#----
    aMethods<-listBuiltInType01()
	mainArgs<-list(clusterFunction="pam", clusterArgs=list(k=3))
	subArgs<-list(resamp.num=20, clusterArgs=list(alpha=0.1))
	for(cf in aMethods){
		set.seed(1045)
		#In sample
	    if(!cf%in%c("hierarchical01","tight")) 
			expect_silent(clusterSingle(biggerMat,  
                subsample=TRUE, 
                sequential=FALSE, 
                makeMissingDiss=TRUE,
                subsampleArgs=c(list(clusterFunction=cf,classifyMethod="InSample"),subArgs), 
                mainClusterArgs=mainArgs,
                isCount=FALSE))
		else expect_message(clusterSingle(biggerMat,  
            subsample=TRUE, 
            sequential=FALSE, 
            makeMissingDiss=TRUE,
            subsampleArgs=c(list(clusterFunction=cf,classifyMethod="InSample"),subArgs), 
            mainClusterArgs=mainArgs),"Making nxn dissimilarity matrix")
       if(!is.null(getBuiltInFunction(cf)@classifyFUN)){
		#All 
   		set.seed(1045)
	    expect_silent(clusterSingle(biggerMat,  subsample=TRUE, sequential=FALSE, subsampleArgs=c(list(clusterFunction=cf,classifyMethod="All"),subArgs), mainClusterArgs=mainArgs,isCount=FALSE))
		#OutOfSample
		set.seed(1045)
	    expect_silent(clusterSingle(biggerMat,  subsample=TRUE, sequential=FALSE, subsampleArgs=c(list(clusterFunction=cf,classifyMethod="OutOfSample"),subArgs), mainClusterArgs=mainArgs,isCount=FALSE))

       }
	}

    ## get NA values
	#set.seed(1045)
	# with the new version, we fix NA's in subsampleClustering
#     expect_error(clusterSingle(mat,
#        subsample=TRUE, sequential=FALSE,
#        subsampleArgs=list(resamp.num=20,clusterArgs=list(k=3),clusterFunction="pam",
# 	   classifyMethod="OutOfSample"),
#        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),isCount=FALSE),"NA values found in dissimilarity matrix")

    # Check warnings in missing args in subsample -- should borrow from mainClusterArgs .
	set.seed(1045)
    expect_warning(clusterSingle(mat,  subsample=TRUE, sequential=FALSE,
		subsampleArgs= list(clusterFunction="pam", resamp.num=3), 
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
		isCount=FALSE),
        "missing arguments k provided from those in 'mainClusterArgs'")
	# Check warnings in missing clusterFunction in subsample -- should borrow from mainClusterArgs .
	expect_message(clusterSingle(mat,  subsample=TRUE, sequential=FALSE,
		subsampleArgs=list(resamp.num=3,clusterArgs=list(k=3)), 
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
		isCount=FALSE),
		"a clusterFunction was not given for subsampleClustering")
	#Check different function types -- should error out.
    expect_error(clusterSingle(mat,  subsample=TRUE, sequential=FALSE, subsampleArgs=list(clusterFunction="pam",resamp.num=3),  mainClusterArgs=list(clusterFunction="tight",clusterArgs=list(alpha=0.1)), isCount=FALSE),
	"must supply arguments: k")


})

test_that("Interactions of mainClustering and subsampling options",{
    ##This once creates unexpected errors even on master, because subsampleArgs not given (even though mainClusterArgs has k, because k isn't the right one for clusterFunction of main (hierarchical01), they aren't passed on to subsample. Unclear if this is the best option...)
    expect_error(ce<-clusterSingle(se, 
         isCount = FALSE, 
         reduceMethod="PCA",
         nDims=10,
         subsample=TRUE,
         sequential=FALSE,
         mainClusterArgs=list( clusterFunction="hierarchical01",minSize=1, clusterArgs=list(alpha=0.1,k = 4))
         )
         ,"must supply arguments: k"
    )
    #check totally weird names passed are caught
    expect_warning(ce<-clusterSingle(se, 
         isCount = FALSE, 
         reduceMethod="PCA",
         nDims=10,
         subsample=TRUE,
         sequential=FALSE,
         subsampleArgs=list(clusterArgs=list(k=3)),
         mainClusterArgs=list( clusterFunction="hierarchical01",myFunkyName=1, clusterArgs=list(alpha=0.1))
         )
         ,"Some arguments passed via mainClusterArgs in mainClustering step do not match the algorithmType of the given ClusterFunction object: myFunkyName"
    )
    # check: if provide diss and mainClust okay with that but 
    # subsampling not, should get error (even if say makeDiss)
	expect_error(clustSubsample <- clusterSingle(dissMat,  
        inputType="diss",
        subsample=TRUE, 
        sequential=FALSE, makeMissingDiss=TRUE, warnings=TRUE,
        subsampleArgs=list(resamp.num=3,
            clusterFunction="kmeans",clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="hierarchical01", 
            clusterArgs=list(alpha=.3)), 
        isCount=FALSE),
        "In subsampling clustering step")
    # check: should be okay even though main doesn't take type X
    # Should NOT trigger calculation of diss matrix.
	expect_silent(clustSubsample <- clusterSingle(mat,  
        inputType="X",
        subsample=TRUE, 
        sequential=FALSE, makeMissingDiss=TRUE, warnings=TRUE,
        subsampleArgs=list(resamp.num=3, clusterFunction="kmeans",clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="hierarchical01", clusterArgs=list(alpha=.3)), 
        isCount=FALSE))
    # check: should be okay but must calculate diss
	expect_message(clustSubsample <- clusterSingle(mat,  
        inputType="X",
        subsample=TRUE, 
        sequential=FALSE, makeMissingDiss=TRUE, warnings=TRUE,
        subsampleArgs=list(resamp.num=3, classifyMethod="InSample",
             clusterFunction="hierarchicalK",clusterArgs=list(k=3)), 
        mainClusterArgs=list(clusterFunction="hierarchical01",
             clusterArgs=list(alpha=.3)), 
        isCount=FALSE),"Making nxn dissimilarity matrix")

})

test_that("Different passed options of seqCluster",{
    #check sequential
    expect_silent(clustSeq <- clusterSingle(mat,
        subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE))
        )
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE), 
        "if sequential=TRUE, must give seqArgs so as to identify k0 and beta")

        
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,mainClusterArgs=list(clusterFunction="pam"),isCount=FALSE,seqArgs=list(k0=5,verbose=FALSE)), "required argument 'beta' is missing")
	expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,mainClusterArgs=list(clusterFunction="pam"),isCount=FALSE,seqArgs=list(beta=0.9,verbose=FALSE)), "required argument 'k0' is missing")

	#right clusterFunctions
	expect_error(clusterSingle(mat, 
        mainClusterArgs=list(clusterFunction="kmeans"), 
        subsample=TRUE, sequential=TRUE, 
        subsampleArgs=list(clusterFunction="pam",n.sample=40), 
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
        "If choosing subsample=TRUE, the clusterFunction used in the mainClustering step must take input that is a categorical")
	expect_error(clusterSingle(dissMat, inputType="diss", 
        mainClusterArgs=list(clusterFunction="tight"),  
        subsample=FALSE, sequential=TRUE, 
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
		"if subsample=FALSE, sequentical clustering can only be implemented with a clusterFunction for mainClustering step with algorithmType 'K'")
	#warning if try to set k
	expect_warning(clusterSingle(mat, 
        mainClusterArgs=list(clusterFunction="pam"), 
        subsample=TRUE, sequential=TRUE, 
        subsampleArgs=list(clusterFunction="pam",
            n.sample=40,clusterArgs=list(k=3)), 
        isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
		"Setting 'k' in subsampleArgs when sequential=TRUE is called will have no effect.")
	expect_warning(clusterSingle(mat, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)), 
        subsample=FALSE, sequential=TRUE, 
        subsampleArgs=list(clusterFunction="pam",n.sample=40), 
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE)),
		"Setting 'k' in mainClusterArgs when sequential clustering is requested will have no effect.")

	#check all algorithms
    kMethods<-listBuiltInTypeK()
	for(cf in kMethods){
		#check if no subsampling
		expect_silent(clusterSingle(dissMat, inputType=inputType(cf)[1],
            mainClusterArgs=list(clusterFunction=cf),
	        subsample=FALSE, sequential=TRUE,
	        isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)))
	}
	kMethods<-listBuiltInType01()
	for(cf in kMethods){
		#check if no subsampling
		expect_silent(clusterSingle(dissMat, inputType=inputType(cf)[1], 
            mainClusterArgs=list(clusterFunction=cf),
	        subsample=TRUE, sequential=TRUE,
			subsampleArgs=list(clusterFunction="pam",n.sample=40),
	        isCount=FALSE,seqArgs=list(k0=5,beta=0.9,verbose=FALSE)))
	}

    ## Check warnings for esoteric arguments
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE,top.can=5.5)), 
        "If top.can parameter of seqArgs is not a whole number, must be between 0 and 1.")
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE,top.can=-1)), 
        "top.can parameter of seqArgs cannot be NA or negative value")
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE,remain.n=3.4)), 
        "remain.n argument of seqArgs must be a positive whole number")
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE,k.min=3.4)), 
        "k.min argument of seqArgs must be a positive whole number")
    expect_error(clusterSingle(mat,subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE,k.max=20.5)), 
        "k.max argument of seqArgs must be a positive whole number")
    expect_silent(clusterSingle(mat,subsample=FALSE, sequential=TRUE,
        mainClusterArgs=list(clusterFunction="pam"),
        isCount=FALSE,
        seqArgs=list(k0=5,beta=0.9,verbose=FALSE,top.can=0.1)))
        
})

test_that("Different direct options of `clusterSingle` ", {
  #check isCount
  expect_silent(clusterSingle(smSimCount,
      subsample=FALSE, sequential=FALSE,
      mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"),
      isCount=TRUE) )
  #suppressWarnings because in addition to error, R prints warning about NAs
  expect_error(suppressWarnings(clusterSingle(smSimData,
                          subsample=FALSE, sequential=FALSE,
                          mainClusterArgs=list(clusterArgs=list(k=3), clusterFunction="pam"),isCount=TRUE)),"User-supplied `transFun` produces NA values",info="test error handling for isCount=TRUE when can't take log")


})

test_that("`clusterSingle` preserves the colData and rowData of SE", {
  expect_silent(cl<-clusterSingle(se,
      subsample=FALSE, sequential=FALSE,
      mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
      isCount=FALSE))

  expect_equal(colData(cl),colData(se))
  expect_equal(rownames(cl),rownames(se))
  expect_equal(colnames(cl),colnames(se))
  expect_equal(metadata(cl),metadata(se))
  expect_equal(rowData(cl),rowData(se))

})

test_that("`clusterSingle` works with parallel subsampling", {
	## Windows does not support mclapply
	skip_on_os("windows")
	expect_silent(clustSubsample2 <- 
		clusterSingle(mat,  subsample=TRUE,sequential=FALSE,
        subsampleArgs=list(clusterFunction="kmeans",resamp.num=3,
        clusterArgs=list(k=3),ncores=2),
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
		isCount=FALSE))

    expect_silent(clustSubsample1 <- 
		clusterSingle(mat,  subsample=TRUE,sequential=FALSE,
        subsampleArgs=list(clusterFunction="kmeans",resamp.num=3,
        clusterArgs=list(k=3),ncores=1),
        mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=3)),
        isCount=FALSE))
})

