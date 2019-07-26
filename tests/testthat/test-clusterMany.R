context("clusterMany")


test_that("`clusterMany` works with matrix, list of data, ClusterExperiment objects, and
          SummarizedExperiments", {
			  #check all builtin methods
    expect_silent(clustNothing <- clusterMany(mat, 
		ks=c(3,4),clusterFunction=listBuiltInFunctions(),
        subsample=FALSE, sequential=FALSE,
        isCount=FALSE,verbose=FALSE, makeMissingDiss=TRUE))
	expect_silent(clustDF <- clusterMany(data.frame(mat),
		ks=c(3,4),clusterFunction=listBuiltInFunctions(),
		subsample=FALSE, sequential=FALSE,
		isCount=FALSE,verbose=FALSE, makeMissingDiss=TRUE))
		   
    expect_is(clustNothing, "ClusterExperiment")
    expect_is(clustNothing, "SingleCellExperiment")
    expect_is(clustNothing, "SummarizedExperiment")

    expect_silent(clustNothing2 <- clusterMany(se, ks=c(3,4),clusterFunction="pam",
       subsample=FALSE, sequential=FALSE,
       isCount=FALSE,verbose=FALSE))
    expect_equal(colData(clustNothing2),colData(se))
    expect_equal(rownames(clustNothing2),rownames(se))
    expect_equal(colnames(clustNothing2),colnames(se))
    expect_equal(metadata(clustNothing2),metadata(se))
    expect_equal(rowData(clustNothing2),rowData(se))

    expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing))
    expect_true(all(clusterTypes(clustNothing)=="clusterMany"))

    #test running on ClusterExperiment Object -- should add the new clustering
    expect_silent(clustNothing3 <- clusterMany(ccSE, 
        ks=c(3,4),clusterFunction="pam",
        subsample=FALSE, sequential=FALSE, verbose=FALSE))
    expect_true(nClusterings(clustNothing3) == nClusterings(ccSE) + 2)
    expect_equal(colData(clustNothing3),colData(ccSE))
    expect_equal(rownames(clustNothing3),rownames(ccSE))
    expect_equal(colnames(clustNothing3),colnames(ccSE))
    expect_equal(metadata(clustNothing3),metadata(ccSE))
    expect_equal(rowData(clustNothing3),rowData(ccSE))
    
    expect_silent(test <- clusterSingle(se,  subsample=FALSE, 
        sequential=FALSE, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=4)),
        isCount=FALSE))
    expect_silent(clustNothing3<- clusterMany(test, 
        ks=c(3,4),clusterFunction="pam",
        subsample=FALSE, sequential=FALSE, verbose=FALSE))
    expect_silent(clustNothing4<- clusterMany(clustNothing3, 
        ks=c(3:4),clusterFunction="pam",
        subsample=FALSE, sequential=FALSE,verbose=FALSE, eraseOld=TRUE))
    expect_equal(clustNothing3,clustNothing4)

    expect_silent(clustNothing5<- clusterMany(clustNothing3, 
        ks=c(5:6),clusterFunction="pam",
        subsample=FALSE, sequential=FALSE,verbose=FALSE, eraseOld=FALSE))
    expect_equal(NCOL(clusterMatrix(clustNothing5)),5)

    expect_silent(ppIndex<-workflowClusterDetails(clustNothing5))
    expect_equal(as.numeric(table(ppIndex[,"iteration"])),c(2,2))
})


test_that("`clusterMany` works with SingleCellExperiment", {
  #check with sce that has dimRed:
  #takes a while with all functions, but sometimes turn up surprises.
  #for some reason if do clusterFunction=listBuiltInFunctions(), 
  #    expect_silent gets warning, but not if run myself.
  for(kk in 1:length(listBuiltInFunctions())){
	  expect_silent(clustNothing2 <- clusterMany(sceSimDataDimRed,
		   ks=c(3,4),clusterFunction=listBuiltInFunctions()[[kk]],
           makeMissingDiss=TRUE,
	       subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))  	
  }
  expect_equal(colData(clustNothing2),colData(sceSimDataDimRed))
  expect_equal(rownames(clustNothing2),rownames(sceSimDataDimRed))
  expect_equal(colnames(clustNothing2),colnames(sceSimDataDimRed))
  expect_equal(metadata(clustNothing2),metadata(sceSimDataDimRed))
  expect_equal(rowData(clustNothing2),rowData(sceSimDataDimRed))
  expect_equal(reducedDims(clustNothing2),reducedDims(sceSimDataDimRed))



})

# # FIXME: need to make this operational
# test_that("`clusterMany` matches `clusterSingle` results", {
#
#     subArgs<-list(clusterFunction="kmeans",classifyMethod="All")
#     expect_message(ce<-clusterMany(smSimSE,
#         isCount = FALSE,
#         ks = 4:5,
#         alphas=c(0.1),
#         reduceMethod="PCA",
#         nReducedDims=10,
#         minSizes=1,
#         subsample=TRUE,
#         sequential=TRUE,
#         clusterFunction="hierarchical01",
#         subsampleArgs=subArgs,
#         makeMissingDiss=TRUE, verbose=FALSE,
#         random.seed=176201),"Not all of the methods requested in 'reduceMethod' have been calculated")
#
#     mainArgs<-list( clusterFunction="hierarchical01",minSize=1,
#         clusterArgs=list(alpha=0.1))
#     seqArgs<-list(k0=4,beta=0.9)
#     set.seed(176201)
#     expect_silent(ce1<-clusterSingle(smSimSE,
#         isCount = FALSE,
#         saveSubsamplingMatrix=TRUE,
#         reduceMethod="PCA",
#         nDims=10,
#         subsample=TRUE,
#         sequential=TRUE,
#         mainClusterArgs=mainArgs,
#         subsampleArgs=subArgs,
#         seqArgs=seqArgs,
#         makeMissingDiss=TRUE
#         ))
#     expect_equal(primaryCluster(ce),primaryCluster(ce1))
#
# })
test_that("`clusterMany` works with reduceMethod a reducedDims",{
  #check picking all dims in single reduceMethod same as apply directly to matrix 
    expect_silent(clustNothing <-
        clusterMany(t(reducedDims(sceSimDataDimRed)[["PCA"]]),
            ks=c(3,4),clusterFunction="pam", 
            reduceMethod="none",
            subsample=FALSE, sequential=FALSE, 
            isCount=FALSE,verbose=FALSE))
    expect_silent(clustNothing3 <- 
        clusterMany(sceSimDataDimRed, 
            ks=c(3,4),clusterFunction="pam", 
            reduceMethod="PCA",
            subsample=FALSE, sequential=FALSE, 
            isCount=FALSE,verbose=FALSE))
    expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))
    expect_equal(NCOL(clusterMatrix(clustNothing)),2)
    expect_equal(NCOL(clusterMatrix(clustNothing3)),2)
    expect_equal(NCOL(clusterMatrix(clustNothing3)),2)
    expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))

    # check picking certain dims in single 
    # reduceMethod same as apply directly to matrix 
    # and that get right reducedDim returned
    expect_silent(clustNothing <- 
        clusterMany(simData, 
            ks=c(3,4),nReducedDims=c(5:6),
            reduceMethod="PCA",
            clusterFunction="pam", subsample=FALSE, 
            sequential=FALSE, isCount=FALSE,verbose=FALSE))
    expect_equal(NCOL(clusterMatrix(clustNothing)),4)
    expect_equal(abs(reducedDim(clustNothing,"PCA")),
        abs(reducedDim(sceSimDataDimRed,"PCA")[,1:6]))

    expect_silent(clustNothing2 <- 
        clusterMany(sceSimDataDimRed, 
        ks=c(3,4),nReducedDims=c(5:6),
        reduceMethod="PCA",
        clusterFunction="pam", subsample=FALSE, 
        sequential=FALSE, isCount=FALSE,verbose=FALSE))
    expect_equal(NCOL(clusterMatrix(clustNothing2)),4)
    expect_equal(reducedDim(clustNothing2,"PCA"), 
        reducedDim(sceSimDataDimRed,"PCA"))

    expect_message(clustNothing4 <- 
        clusterMany(sceSimData, 
            ks=c(3,4),nReducedDims=c(5:6),
            reduceMethod="PCA", clusterFunction="pam", 
            subsample=FALSE, sequential=FALSE, 
            isCount=FALSE,verbose=FALSE),
            "Not all of the methods requested in 'reduceMethod' have been calculated")
    expect_equal(clusterMatrix(clustNothing2), 
        clusterMatrix(clustNothing4))
  
    expect_message(clustNothingBoth1<-
        clusterMany(clustNothing4,
            ks=c(3,4),nReducedDims=c(5:6),
            reduceMethod=c("PCA","mad"),
            clusterFunction="pam", subsample=FALSE, 
            sequential=FALSE, verbose=FALSE),
            "Not all of the methods requested in 'reduceMethod' have been calculated")
    expect_silent(clustNothingBoth2<-
        clusterMany(clustNothingBoth1, 
            ks=c(3,4),nReducedDims=c(5:6),
            reduceMethod=c("PCA","mad"),
            clusterFunction="pam", subsample=FALSE, 
            sequential=FALSE, verbose=FALSE))
    expect_equal(clusterMatrix(clustNothingBoth1,
        which="clusterMany"), clusterMatrix(clustNothingBoth2,
        which="clusterMany"))

    #check picking reduceMethod="none" same as apply directly to matrix 
    expect_silent(clustNothing <- 
        clusterMany(simData, ks=c(3,4),
            clusterFunction="pam", subsample=FALSE, 
            sequential=FALSE, isCount=FALSE,verbose=FALSE))
    expect_silent(clustNothing3 <- 
        clusterMany(sceSimDataDimRed, ks=c(3,4),
            clusterFunction="pam", reduceMethod="none",
            subsample=FALSE, sequential=FALSE, 
            isCount=FALSE,verbose=FALSE))
    expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))

    #checks that nReducedDims ignored if reduceMethod="none"
    expect_silent(clustNothing <- 
        clusterMany(simData,
            ks=c(3,4),nReducedDims=c(5:6),
            reduceMethod="none",
            clusterFunction="pam", subsample=FALSE, 
            sequential=FALSE, isCount=FALSE,
            verbose=FALSE))
    expect_silent(clustNothing3 <- 
        clusterMany(sceSimDataDimRed,
            ks=c(3,4),nReducedDims=c(5:6),
            reduceMethod="none",
            clusterFunction="pam", subsample=FALSE, 
            sequential=FALSE, isCount=FALSE,verbose=FALSE))
    expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))

  
})
test_that("`clusterMany` works with reduceMethod a filtering",{
  expect_message(clustNothing4 <- clusterMany(sceSimDataDimRed, 
                                             ks=c(3,4),clusterFunction="pam", reduceMethod="mad",
                                             subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE),"Not all of the methods requested in 'reduceMethod' have been calculated")
  
  expect_error(clustNothing5 <- clusterMany(sceSimDataDimRed, 
                                             ks=c(3,4),clusterFunction="pam", reduceMethod=c("mad","cv"),
                                             subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE),"do not match any stored or built-in filtering statistics or dimensionality reduction")
  expect_message(clustNothing5 <- clusterMany(sceSimDataDimRed, 
                                            ks=c(3,4),clusterFunction="pam", reduceMethod=c("mad","abscv"),
                                            subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE),"Not all of the methods requested in 'reduceMethod' have been calculated")
  expect_true(all(filterNames(sceSimDataDimRed) %in% filterNames(clustNothing5)))
  expect_true(all(c("mad","abscv","var","mean") %in% filterNames(clustNothing5)))
  
  expect_silent(clusterMany(clustNothing5, 
                                              ks=c(3,4),clusterFunction="pam", reduceMethod=c("mad","abscv"),
                                              subsample=FALSE, sequential=FALSE, verbose=FALSE))
})
test_that("`clusterMany` works with hdf5", {
	########
	#Check if use PCA (not in hdf5) changes nothing, as expect.
	########
    for(kk in 1:length(listBuiltInFunctions())){
  	  expect_silent(clustNothing2 <- clusterMany(hdfSCE,
  		   ks=c(3,4),clusterFunction=listBuiltInFunctions()[[kk]],
  	       makeMissingDiss=TRUE,
           subsample=FALSE, sequential=FALSE, 
           isCount=FALSE,verbose=FALSE))  	
    }
    expect_silent(clustNothing <- 
        clusterMany(t(reducedDims(sceSimDataDimRed)[["PCA"]]), 
    	    ks=c(3,4),clusterFunction="pam", 
            reduceMethod="none",
            subsample=FALSE, sequential=FALSE, 
            isCount=FALSE,verbose=FALSE))
    expect_silent(clustNothing3 <- 
        clusterMany(hdfSCE, 
            ks=c(3,4),clusterFunction="pam", 
            reduceMethod="PCA", subsample=FALSE, 
            sequential=FALSE, isCount=FALSE,verbose=FALSE))
    expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))
    expect_equal(NCOL(clusterMatrix(clustNothing)),2)
    expect_equal(NCOL(clusterMatrix(clustNothing3)),2)

	########
	#Check if not using PCA (reduceMethod = "none")
	########
    for(kk in 1:length(listBuiltInFunctions())){
  	  expect_silent(clustNothing2 <- 
          clusterMany(hdfSCE,reduceMethod = "none", 
  		   ks=c(3,4), makeMissingDiss=TRUE,
           clusterFunction=listBuiltInFunctions()[[kk]],
  	       subsample=FALSE, sequential=FALSE, 
           isCount=FALSE,verbose=FALSE))  	
    }

	########
	#Check directly on hdf5 object
	########
    for(kk in 1:length(listBuiltInFunctions())){
  	  expect_silent(clustNothing2 <- 
          clusterMany(assay(hdfObj),reduceMethod = "none", 
  		       ks=c(3,4),clusterFunction=listBuiltInFunctions()[[kk]],
  	           subsample=FALSE, sequential=FALSE, 
               makeMissingDiss=TRUE,
               isCount=FALSE,verbose=FALSE))  	
    }

})

test_that("`clusterMany` works changing parameters", {
	#--------
  #check dim reduce in combination
	#--------
    expect_silent(cc <- clusterMany(mat, 
        ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),
  	    reduceMethod=c("none","PCA","var","abscv","mad"),
        clusterFunction="pam",
        subsample=FALSE, sequential=FALSE,
        verbose=FALSE, isCount=FALSE))
	expect_equal(sort(reducedDimNames(cc)),sort(c("PCA")))
	expect_equal(sort(filterNames(cc)),sort(c("var","abscv","mean","mad")))

	expect_message(cc2 <- clusterMany(se, 
        ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),
	    reduceMethod=c("none","PCA","var","abscv","mad"),
        clusterFunction="pam",
        subsample=FALSE, sequential=FALSE,
        verbose=FALSE, isCount=FALSE),
        "Not all of the methods requested in 'reduceMethod' have been calculated")
    expect_equal(sort(reducedDimNames(cc2)),sort(c("PCA")))
    expect_equal(sort(filterNames(cc2)),
        sort(c("b","var","abscv","mean","mad")))

	expect_message(cc3 <- clusterMany(sce, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),
	reduceMethod=c("none","PCA","var","abscv","mad"),clusterFunction="pam",
  subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),"Not all of the methods requested in 'reduceMethod' have been calculated")
expect_equal(sort(reducedDimNames(cc3)),sort(c("PCA")))
expect_equal(sort(filterNames(cc3)),sort(c("b","var","abscv","mean","mad")))

	#Only existing values 
	expect_silent(cc4 <- clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
	reduceMethod=c("none","Red1","Filter1","Filter2"),clusterFunction="pam",
  subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
expect_equal(sort(reducedDimNames(cc4)),sort(c("Red1")))
expect_equal(sort(filterNames(cc4)),sort(c("b","Filter1","Filter2")))

  expect_silent(ceReRun <- clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
  	reduceMethod=c("none","Red1","Filter1","Filter2"),clusterFunction="pam",
    subsample=FALSE, sequential=FALSE,verbose=FALSE))
	expect_equal(sort(reducedDimNames(ceReRun)),sort(reducedDimNames(cc4)))
	expect_equal(sort(filterNames(ceReRun)),sort(filterNames(cc4)))


	#Only existing values 
	expect_silent(cc4 <- clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
	reduceMethod=c("none","Red1","Filter1","Filter2"),clusterFunction="pam",
  subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
expect_equal(sort(reducedDimNames(cc4)),sort(c("Red1")))
expect_equal(sort(filterNames(cc4)),sort(c("b","Filter1","Filter2")))


	#--------
	# Mixing saved and unsaved (gives warnings/errors)
	#--------
	#following gives warning because can't mix saved and calculate internally
	notAllWarning<-"All values of 'reduceMethod' need to either match an existing"
	expect_error(clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
			reduceMethod=c("PCA","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),notAllWarning)
	expect_error(clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
			reduceMethod=c("var","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),notAllWarning)
	expect_error(clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
			reduceMethod=c("PCA","Filter1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),notAllWarning)
	expect_error(clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("var","Filter1"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),
		notAllWarning)	
	#repeat for ce 
	expect_error(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
			reduceMethod=c("PCA","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE),notAllWarning)
	expect_error(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
			reduceMethod=c("PCA","Filter1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE),notAllWarning)
	expect_error(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
			reduceMethod=c("var","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE),notAllWarning)
	expect_error(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("var","Filter1"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE),
		notAllWarning)	
	
	#--------
	# Mixing mix across filter/reduceMethod okay
	#--------
	sceSimData2<-sceFull
	reducedDims(sceSimData2)<-SimpleList()
	expect_message(c1<-clusterMany(sceSimData2, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("PCA"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),"Not all of the methods requested in 'reduceMethod' have been calculated.")
	expect_equal(clusterExperiment:::filterStats(c1),clusterExperiment:::filterStats(sceSimData2))
	expect_equal(reducedDimNames(c1),"PCA")
	
	sceSimData3<-sceFull
	rowData(sceSimData3)<-NULL
	expect_message(c2<-clusterMany(sceSimData3, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("var"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),"Not all of the methods requested in 'reduceMethod' have been calculated.")
	expect_equal(reducedDims(c2),reducedDims(sceSimData3))
	expect_equal(filterNames(c2),"var")
	
	#repeat for ce
	ce5<-cc4
	reducedDims(ce5)<-SimpleList()
	expect_message(c1<-clusterMany(ce5, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("PCA"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE),"Not all of the methods requested in 'reduceMethod' have been calculated.")
	expect_equal(clusterExperiment:::filterStats(c1),clusterExperiment:::filterStats(ce5))
	expect_equal(reducedDimNames(c1),"PCA")
	
	ce6<-cc4
	rowData(ce6)<-NULL
	expect_message(c2<-clusterMany(ce6, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("var"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE),"Not all of the methods requested in 'reduceMethod' have been calculated.")
	expect_equal(reducedDims(c2),reducedDims(ce6))
	expect_equal(filterNames(c2),"var")

	
	#--------
	#check if only dim reduction
	#--------
	expect_silent(clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("Red1"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	sceFull2<-makeReducedDims(sceFull,reducedDims="PCA",maxDims=10,isCount=FALSE)
	expect_silent(clusterMany(sceFull2, ks=c(3,4),nFilterDims=c(10),nReducedDims=c(2),
			reduceMethod=c("PCA","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	expect_silent(clusterMany(sceFull2, ks=c(3,4),nFilterDims=c(10),nReducedDims=c(2),
			reduceMethod=c("Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))

	#--------
	#check only filter
	#--------
	expect_silent(clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("Filter1"),clusterFunction="pam",
  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	sceFull2<-makeFilterStats(sceFull,filterStat="var",isCount=FALSE)
	
	expect_silent(clusterMany(sceFull2, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("var","Filter1"),clusterFunction="pam",
  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	expect_silent(clusterMany(sceFull2, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("var"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	
	#--------
	#check only none
	#--------
	expect_silent(clusterMany(sceFull, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(2),
		reduceMethod=c("none"),clusterFunction="pam",
  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))

	#--------
	#check returning paramMatrix
	#--------
	expect_silent(param <- clusterMany(mat, 
	  ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),reduceMethod=c("none","PCA","var"),
	  clusterFunction="pam", subsample=FALSE, sequential=FALSE,run=FALSE,verbose=FALSE,
	  isCount=FALSE))
  # #check that giving param$paramMatrix works -- not implemented...
  # cc2 <- clusterMany(mat, ks=c(3,4),nFilterDims=c(10, 15),nReducedDims=c(3,4),
  # 	reduceMethod=c("none","PCA","var"),clusterFunction="pam",
  #   subsample=FALSE, sequential=FALSE,verbose=FALSE,
  #    isCount=FALSE,paramMatrix=param$paramMatrix, mainClusterArgs=param$mainClusterArgs,
  # 	 seqArgs=param$seqArgs,subsampleArgs=param$subsampleArgs)
  # expect_equal(cc,cc2)
  
  expect_silent(cc2 <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
                    distFunction=c("euclidean","manhattan",NA),
                    subsample=FALSE, sequential=FALSE,verbose=FALSE,
                    isCount=FALSE, makeMissingDiss=TRUE))

  # # #check giving distance -- this still doesn't work.
  # # #problem to just give dist because need to grab from global environment
  # # # So doesn't work in testthat
  # dist1<-function(x){dist(x,method="manhattan")}
  # dist2<-function(x){dist(x)}
  # expect_silent(cc <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
  #                   distFunction=c("dist1","dist2",NA),
  #                   subsample=FALSE, sequential=FALSE,verbose=FALSE,
  #                   isCount=FALSE))
  # expect_equal(unname(clusterMatrix(cc)),unname(clusterMatrix(cc2)))
  # 


  #check doesn't spit out warnings because alphas/mainClustering args not match 
  expect_silent(clusterMany(mat, 
      clusterFunction=c("pam","hierarchical01"),
      ks=c(3,4),
      alphas=c(0.1,0.2),
      subsample=FALSE, sequential=FALSE,
      verbose=FALSE,
      mainClusterArgs=list(clusterArgs=list(evalClusterMethod="average")),
      isCount=FALSE))
  
  #check doesn't spit out warnings because alphas/mainClustering args not match 
  expect_silent(clusterMany(mat, 
      clusterFunction=c("pam","hierarchical01"),ks=c(3,4),
      betas=c(.7,.9), minSizes=c(3,5),
      subsample=FALSE, sequential=FALSE,verbose=FALSE,
      mainClusterArgs=list(clusterArgs=list(evalClusterMethod="average")),
      isCount=FALSE))
})

test_that("`getClusterManyParams` works", {
	cc<-clusterMany(mat, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),
		reduceMethod=c("none","PCA","var"),clusterFunction="pam",
	  	subsample=FALSE, sequential=FALSE,run=TRUE,verbose=FALSE,
	    isCount=FALSE)
	cc<-makeConsensus(cc,proportion=1,whichClusters = "clusterMany")
	expect_silent(paramAll<-getClusterManyParams(cc))
	expect_equal(colnames(paramAll),c("clusteringIndex", "reduceMethod", "nReducedDims", "nFilterDims", "k"))
	expect_true(is.data.frame(paramAll))
	expect_equal(sort(as.character(unique(paramAll[,"reduceMethod"]))),sort(c("none","var","PCA")))
	
	expect_equal(sort(unique(paramAll[,"nReducedDims"]),na.last=TRUE),sort(c(NA,3,4),na.last=TRUE))
	expect_equal(sort(unique(paramAll[,"nFilterDims"]),na.last=TRUE),sort(c(NA,10,15),na.last=TRUE))
	expect_true(is.numeric(paramAll[,"k"]))
	
    # check simplify=TRUE
	expect_silent(paramSub<-getClusterManyParams(cc,whichClusters=c(4,6),
        simplify=TRUE))
	expect_equal(colnames(paramSub),c("clusteringIndex", "reduceMethod","nFilterDims"))
	expect_true(is.data.frame(paramSub))
	expect_equal(sort(unique(paramSub[,"nFilterDims"]),na.last=TRUE),c(10,NA))
	
	expect_warning(getClusterManyParams(cc,whichClusters=1:5),"some clusters indicated in 'whichClusters' do not have type 'clusterMany'")
	expect_warning(getClusterManyParams(cc,whichClusters=1),"did not return any clusters of type 'clusterMany'")
	expect_warning(getClusterManyParams(cc,whichClusters="mergeClusters"),"did not return any clusters")
})



test_that("`clusterMany` consistent results (no transformation)", {
	
	#----
	#check get the same result from previous times when transform the data
	#Relies on having these checks in clusterSingle, so know clusterSingle is right
    #----
	#make it big enough can do pca and filter...
	contData<-simData[,1:20]
	expectTrans1<-round(contData[1,],2)
    testSE<-SummarizedExperiment(contData)
    testSCE<-as(testSE,"SingleCellExperiment")
    
	#matrix
	expect_silent(ccVar<-clusterSingle(contData,
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
        isCount=FALSE))
   	expect_silent(ccPCA<-clusterSingle(contData, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
        isCount=FALSE))
  	expect_silent(ccNone<-clusterSingle(contData, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
        isCount=FALSE))
	expect_silent(cm<-clusterMany(contData, clusterFunction="pam",ks=3,
        subsample=FALSE, sequential=FALSE, 
        reduceMethod=c("PCA","var","none"), verbose=FALSE,
        nReducedDims=3, nFilterDims=3,isCount=FALSE))
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
    
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 
	

   #SE
	expect_silent(ccVar<-clusterSingle(testSE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
 	   isCount=FALSE))
  	expect_silent(ccPCA<-clusterSingle(testSE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
 	   isCount=FALSE))
 	expect_silent(ccNone<-clusterSingle(testSE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
 	   isCount=FALSE))
	expect_message(cm<-clusterMany(testSE, clusterFunction="pam",ks=3,
   	        subsample=FALSE, sequential=FALSE, 
            reduceMethod=c("PCA","var","none"),
   	        nReducedDims=3, nFilterDims=3,isCount=FALSE),
            "Not all of the methods requested in 'reduceMethod' have been calculated.")
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 

	
   #SCE
	expect_silent(ccVar<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
 	   isCount=FALSE))
  	expect_silent(ccPCA<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
 	   isCount=FALSE))
 	expect_silent(ccNone<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA, 
        mainClusterArgs=list(clusterFunction="pam",
            clusterArgs=list(k=3)),
 	   isCount=FALSE))
	expect_message(cm<-clusterMany(testSCE, clusterFunction="pam",ks=3,
        subsample=FALSE, sequential=FALSE, 
        reduceMethod=c("PCA","var","none"),
        nReducedDims=3, nFilterDims=3,isCount=FALSE),
        "Not all of the methods requested in 'reduceMethod' have been calculated.")
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 
	

   #CE
	expect_silent(ccVar2<-clusterSingle(ccVar, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, clusterLabel="redo",
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3))))
	expect_silent(ccPCA2<-clusterSingle(ccPCA, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, clusterLabel="redo",
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3))))
	expect_silent(ccNone2<-clusterSingle(ccNone, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA,clusterLabel="redo",
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3))))
	expect_silent(cm2<-clusterMany(cm, clusterFunction="pam",ks=3,
        subsample=FALSE, sequential=FALSE, 
        reduceMethod=c("PCA","var","none"), verbose=FALSE,
        nReducedDims=3, nFilterDims=3))
	expect_equal(nClusterings(cm2),6)	
	expect_silent(params<-getClusterManyParams(cm2))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm2)[1,],2), expectTrans1) 
    

  
})


test_that("`clusterMany` consistent results (with transformation)", {
	
	#----
	#check get the same result from previous times when transform the data
	#Relies on tests in clusterSingle being right
    #----
	#make it big enough can do pca and filter...
	countData<-simCount[,1:20]
    testSE<-SummarizedExperiment(countData)
    testSCE<-as(testSE,"SingleCellExperiment")
	expectTrans1<-round(log2(countData[1,]+1),2)

	#matrix
	expect_silent(ccVar<-clusterSingle(countData, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
   	expect_silent(ccPCA<-clusterSingle(countData, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
  	expect_silent(ccNone<-clusterSingle(countData, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
	expect_silent(cm<-clusterMany(countData, clusterFunction="pam",ks=3,
        subsample=FALSE, sequential=FALSE, 
        reduceMethod=c("PCA","var","none"), verbose=FALSE,
        nReducedDims=3, nFilterDims=3,isCount=TRUE))
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
    
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 
	

   #SE
	expect_silent(ccVar<-clusterSingle(testSE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
  	expect_silent(ccPCA<-clusterSingle(testSE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
 	expect_silent(ccNone<-clusterSingle(testSE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
        isCount=TRUE))
	expect_message(cm<-clusterMany(testSE, clusterFunction="pam",ks=3,
        subsample=FALSE, sequential=FALSE, 
        reduceMethod=c("PCA","var","none"),
        nReducedDims=3, nFilterDims=3,isCount=TRUE),"Not all of the methods requested in 'reduceMethod' have been calculated")
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 

	
   #SCE
	expect_silent(ccVar<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
        isCount=TRUE))
  	expect_silent(ccPCA<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, 
        mainClusterArgs=list(clusterFunction="pam",
        clusterArgs=list(k=3)),
 	   isCount=TRUE))
 	expect_silent(ccNone<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA, 
        mainClusterArgs=list(clusterFunction="pam",
        clusterArgs=list(k=3)),
 	   isCount=TRUE))
	expect_message(cm<-clusterMany(testSCE, clusterFunction="pam",ks=3,
        subsample=FALSE, sequential=FALSE,  
        reduceMethod=c("PCA","var","none"),
        nReducedDims=3, nFilterDims=3,isCount=TRUE),"Not all of the methods requested in 'reduceMethod' have been calculated")
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 
	

   #CE
	expect_silent(ccVar2<-clusterSingle(ccVar, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, clusterLabel="redo",
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3))))
	expect_silent(ccPCA2<-clusterSingle(ccPCA, 
        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
        nDims=3, clusterLabel="redo",
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3))))
	expect_silent(ccNone2<-clusterSingle(ccNone, 
        subsample=FALSE, sequential=FALSE, reduceMethod="none",
        nDims=NA,clusterLabel="redo",
		mainClusterArgs=list(clusterFunction="pam",
             clusterArgs=list(k=3))))
	expect_silent(cm2<-clusterMany(cm, clusterFunction="pam",ks=3,
        subsample=FALSE, sequential=FALSE, 
        reduceMethod=c("PCA","var","none"), verbose=FALSE,
        nReducedDims=3, nFilterDims=3))
	expect_equal(nClusterings(cm2),6)	
	expect_silent(params<-getClusterManyParams(cm2))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm2)[1,],2), expectTrans1) 
    
  
})

test_that("`clusterMany` works with ClusterFunction objects",{
		expect_silent(clustAll1<-clusterMany(sceSimDataDimRed,reduceMethod="PCA",
			   ks=c(3,4),clusterFunction=listBuiltInFunctions()[1:2],
		       subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE,random.seed=1250))  	
	  
	expect_silent(clustAll2<-clusterMany(sceSimDataDimRed,reduceMethod="PCA",
	 ks=c(3,4),clusterFunction=getBuiltInFunction(listBuiltInFunctions()[1:2]),
		       subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE,random.seed=1250))  	
	  
		expect_equal(clustAll2,clustAll1)			 
})
