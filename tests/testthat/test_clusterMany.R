context("clusterMany")
source("create_objects.R")

test_that("`clusterMany` works with matrix, list of data, ClusterExperiment objects, and
          SummarizedExperiments", {
			  #check all builtin methods
    expect_silent(clustNothing <- clusterMany(mat, 
		ks=c(3,4),clusterFunction=listBuiltInFunctions(),
        subsample=FALSE, sequential=FALSE,
        isCount=FALSE,verbose=FALSE))
	expect_silent(clustDF <- clusterMany(data.frame(mat),
		ks=c(3,4),clusterFunction=listBuiltInFunctions(),
		subsample=FALSE, sequential=FALSE,
		isCount=FALSE,verbose=FALSE))
		   
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
    expect_silent(clustNothing3 <- clusterMany(ccSE, ks=c(3,4),clusterFunction="pam",
                                 subsample=FALSE, sequential=FALSE,verbose=FALSE))
    expect_true(nClusterings(clustNothing3) == nClusterings(ccSE) + 2)
    expect_equal(colData(clustNothing3),colData(ccSE))
    expect_equal(rownames(clustNothing3),rownames(ccSE))
    expect_equal(colnames(clustNothing3),colnames(ccSE))
    expect_equal(metadata(clustNothing3),metadata(ccSE))
    expect_equal(rowData(clustNothing3),rowData(ccSE))
    
    expect_silent(test <- clusterSingle(se,  subsample=FALSE, sequential=FALSE, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=4)),isCount=FALSE))
    expect_silent(clustNothing3<- clusterMany(test, ks=c(3,4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,verbose=FALSE))
    expect_silent(clustNothing4<- clusterMany(clustNothing3, ks=c(3:4),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,verbose=FALSE, eraseOld=TRUE))
    expect_equal(clustNothing3,clustNothing4)

    clustNothing5<- clusterMany(clustNothing3, ks=c(5:6),clusterFunction="pam",
                                   subsample=FALSE, sequential=FALSE,verbose=FALSE, eraseOld=FALSE)
    expect_equal(NCOL(clusterMatrix(clustNothing5)),5)

    ppIndex<-workflowClusterDetails(clustNothing5)
    expect_equal(as.numeric(table(ppIndex[,"iteration"])),c(2,2))
})
test_that("`clusterMany` works with SingleCellExperiment", {
  #check with sce that has dimRed:
  #takes a while with all functions, but sometimes turn up surprises.
  #for some reason if do clusterFunction=listBuiltInFunctions(), 
  #    expect_silent gets warning, but not if run myself.
  for(kk in 1:length(listBuiltInFunctions)){
	  expect_silent(clustNothing2 <- clusterMany(sceSimDataDimRed,
		   ks=c(3,4),clusterFunction=listBuiltInFunctions()[[kk]],
	       subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))  	
  }
  expect_equal(colData(clustNothing2),colData(sceSimDataDimRed))
  expect_equal(rownames(clustNothing2),rownames(sceSimDataDimRed))
  expect_equal(colnames(clustNothing2),colnames(sceSimDataDimRed))
  expect_equal(metadata(clustNothing2),metadata(sceSimDataDimRed))
  expect_equal(rowData(clustNothing2),rowData(sceSimDataDimRed))
  expect_equal(reducedDims(clustNothing2),reducedDims(sceSimDataDimRed))

  #check picking all dims in single reduceMethod same as apply directly to matrix 
  expect_silent(clustNothing <- clusterMany(t(reducedDims(sceSimDataDimRed)[["PCA"]]), 
  	ks=c(3,4),clusterFunction="pam", reduceMethod="none",
    subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_silent(clustNothing3 <- clusterMany(sceSimDataDimRed, 
	ks=c(3,4),clusterFunction="pam", reduceMethod="PCA",
	subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))
  expect_equal(NCOL(clusterMatrix(clustNothing)),2)
  expect_equal(NCOL(clusterMatrix(clustNothing3)),2)

  #check picking certain dims in single reduceMethod same as apply directly to matrix 
  #and that get right reducedDim returned
  expect_silent(clustNothing <- clusterMany(simData, 
	  ks=c(3,4),nReducedDims=c(5:6),reduceMethod="PCA",
  	  clusterFunction="pam", subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_equal(NCOL(clusterMatrix(clustNothing)),4)
  expect_equal(abs(reducedDim(clustNothing,"PCA")), abs(reducedDim(sceSimDataDimRed,"PCA")[,1:6]))

  expect_silent(clustNothing3 <- clusterMany(sceSimDataDimRed, 
	  ks=c(3,4),nReducedDims=c(5:6),reduceMethod="PCA",
      clusterFunction="pam", subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_equal(NCOL(clusterMatrix(clustNothing3)),4)
  expect_equal(reducedDim(clustNothing3,"PCA"), reducedDim(sceSimDataDimRed,"PCA"))

  expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))

  #check picking reduceMethod="none" same as apply directly to matrix 
  expect_silent(clustNothing <- clusterMany(simData, ks=c(3,4),clusterFunction="pam", subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_silent(clustNothing3 <- clusterMany(sceSimDataDimRed, ks=c(3,4),clusterFunction="pam", reduceMethod="none",subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))

  #checks that nReducedDims ignored if reduceMethod="none"
  expect_silent(clustNothing <- clusterMany(simData, 
	  ks=c(3,4),nReducedDims=c(5:6),reduceMethod="none",
  	  clusterFunction="pam", subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_silent(clustNothing3 <- clusterMany(sceSimDataDimRed, 
	  ks=c(3,4),nReducedDims=c(5:6),reduceMethod="none",
      clusterFunction="pam", subsample=FALSE, sequential=FALSE, isCount=FALSE,verbose=FALSE))
  expect_equal(clusterMatrix(clustNothing), clusterMatrix(clustNothing3))



})

test_that("`clusterMany` works changing parameters", {
	#--------
  #check dim reduce in combination
	#--------
  expect_silent(cc <- clusterMany(mat, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(3,4),
  	reduceMethod=c("none","PCA","var","abscv","mad"),clusterFunction="pam",
    subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	expect_equal(sort(reducedDimNames(cc)),sort(c("PCA")))
	expect_equal(sort(filterNames(cc)),sort(c("var","abscv","mean","mad")))

	expect_silent(cc2 <- clusterMany(se, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(3,4),
	reduceMethod=c("none","PCA","var","abscv","mad"),clusterFunction="pam",
  subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
expect_equal(sort(reducedDimNames(cc2)),sort(c("PCA")))
expect_equal(sort(filterNames(cc2)),sort(c("var","abscv","mean","mad")))

	expect_silent(cc3 <- clusterMany(scf, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(3,4),
	reduceMethod=c("none","PCA","var","abscv","mad"),clusterFunction="pam",
  subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
expect_equal(sort(reducedDimNames(cc3)),sort(c("PCA")))
expect_equal(sort(filterNames(cc3)),sort(c("var","abscv","mean","mad")))

	#Only existing values 
	expect_silent(cc4 <- clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
	reduceMethod=c("none","Red1","Filter1","Filter2"),clusterFunction="pam",
  subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
expect_equal(sort(reducedDimNames(cc4)),sort(c("Red1")))
expect_equal(sort(filterNames(cc4)),sort(c("Filter1","Filter2")))

  expect_silent(ceReRun <- clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
  	reduceMethod=c("none","Red1","Filter1","Filter2"),clusterFunction="pam",
    subsample=FALSE, sequential=FALSE,verbose=FALSE))
	expect_equal(sort(reducedDimNames(ceReRun)),sort(reducedDimNames(cc4)))
	expect_equal(sort(filterNames(ceReRun)),sort(filterNames(cc4)))


	#Only existing values 
	expect_silent(cc4 <- clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
	reduceMethod=c("none","Red1","Filter1","Filter2"),clusterFunction="pam",
  subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
expect_equal(sort(reducedDimNames(cc4)),sort(c("Red1")))
expect_equal(sort(filterNames(cc4)),sort(c("Filter1","Filter2")))


	#--------
	# Mixing saved and unsaved (gives warnings/errors)
	#--------
	#following gives warning because can't mix saved and calculate internally
	expect_warning(clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
			reduceMethod=c("PCA","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),"Not all of reduceMethod value match a reducedDimNames or filterNames")
	expect_warning(clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
			reduceMethod=c("var","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),"Not all of reduceMethod value match a reducedDimNames or filterNames")
	expect_warning(clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
			reduceMethod=c("PCA","Filter1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),"Not all of reduceMethod value match a reducedDimNames or filterNames")
	expect_warning(clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("var","Filter1"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE),
		"Not all of reduceMethod value match a reducedDimNames or filterNames")	
	#repeat for ce 
	expect_warning(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
			reduceMethod=c("PCA","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE),"Not all of reduceMethod value match a reducedDimNames or filterNames")
	expect_warning(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
			reduceMethod=c("PCA","Filter1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE),"Not all of reduceMethod value match a reducedDimNames or filterNames")
	expect_warning(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
			reduceMethod=c("var","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE),"Not all of reduceMethod value match a reducedDimNames or filterNames")
	expect_warning(clusterMany(cc4, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("var","Filter1"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE),
		"Not all of reduceMethod value match a reducedDimNames or filterNames")	
	
	#--------
	# Mixing mix across filter/reduceMethod okay
	#--------
	scf2<-scfFull
	reducedDims(scf2)<-SimpleList()
	expect_silent(c1<-clusterMany(scf2, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("PCA"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	expect_equal(filterStats(c1),filterStats(scf2))
	expect_equal(reducedDimNames(c1),"PCA")
	
	scf3<-scfFull
	filterStats(scf3)<-NULL
	expect_silent(c2<-clusterMany(scf3, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("var"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	expect_equal(reducedDims(c2),reducedDims(scf3))
	expect_equal(filterNames(c2),"var")
	
	#repeat for ce
	ce5<-cc4
	reducedDims(ce5)<-SimpleList()
	expect_silent(c1<-clusterMany(ce5, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("PCA"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE))
	expect_equal(filterStats(c1),filterStats(ce5))
	expect_equal(reducedDimNames(c1),"PCA")
	
	ce6<-cc4
	filterStats(ce6)<-NULL
	expect_silent(c2<-clusterMany(ce6, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("var"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE))
	expect_equal(reducedDims(c2),reducedDims(ce6))
	expect_equal(filterNames(c2),"var")

	
	#--------
	#check if only dim reduction
	#--------
	expect_silent(clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("Red1"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	scfFull2<-makeReducedDims(scfFull,reducedDims="PCA",maxDims=10,isCount=FALSE)
	expect_silent(clusterMany(scfFull2, ks=c(3,4),nFilterDims=c(10),nPCADim=c(2),
			reduceMethod=c("PCA","Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	expect_silent(clusterMany(scfFull2, ks=c(3,4),nFilterDims=c(10),nPCADim=c(2),
			reduceMethod=c("Red1"),clusterFunction="pam",
	  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))

	#--------
	#check only filter
	#--------
	expect_silent(clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("Filter1"),clusterFunction="pam",
  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	scfFull2<-makeFilterStats(scfFull,filterStat="var",isCount=FALSE)
	
	expect_silent(clusterMany(scfFull2, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("var","Filter1"),clusterFunction="pam",
  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	expect_silent(clusterMany(scfFull2, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("var"),clusterFunction="pam",
		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))
	
	#--------
	#check only none
	#--------
	expect_silent(clusterMany(scfFull, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(2),
		reduceMethod=c("none"),clusterFunction="pam",
  		subsample=FALSE, sequential=FALSE,verbose=FALSE, isCount=FALSE))

	#--------
	#check returning paramMatrix
	#--------
	expect_silent(param <- clusterMany(mat, 
	  ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(3,4),reduceMethod=c("none","PCA","var"),
	  clusterFunction="pam", subsample=FALSE, sequential=FALSE,run=FALSE,verbose=FALSE,
	  isCount=FALSE))
  # #check that giving param$paramMatrix works -- not implemented...
  # cc2 <- clusterMany(mat, ks=c(3,4),nFilterDims=c(10, 15),nPCADim=c(3,4),
  # 	reduceMethod=c("none","PCA","var"),clusterFunction="pam",
  #   subsample=FALSE, sequential=FALSE,verbose=FALSE,
  #    isCount=FALSE,paramMatrix=param$paramMatrix, mainClusterArgs=param$mainClusterArgs,
  # 	 seqArgs=param$seqArgs,subsampleArgs=param$subsampleArgs)
  # expect_equal(cc,cc2)
  
  # #check giving distance -- this still doesn't work.
  # #problem to just give dist because need to grab from global environment
  # dist1<-function(x){dist(x,method="manhattan")}
  # dist2<-function(x){dist(x)} 
  # cc <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
  #                   distFunction=c("dist1","dist2",NA),
  #                   subsample=FALSE, sequential=FALSE,verbose=FALSE,
  #                   isCount=FALSE)
  
  #check doesn't spit out warnings because alphas/mainClustering args not match 
  expect_silent(clusterMany(mat, clusterFunction=c("pam","hierarchical01"),ks=c(3,4),
                    alphas=c(0.1,0.2),
                    subsample=FALSE, sequential=FALSE,verbose=FALSE,
                    mainClusterArgs=list(clusterArgs=list(evalClusterMethod="average")),
                    isCount=FALSE))
  
  #check doesn't spit out warnings because alphas/mainClustering args not match 
  expect_silent(clusterMany(mat, clusterFunction=c("pam","hierarchical01"),ks=c(3,4),
                            betas=c(.7,.9), minSizes=c(3,5),
                            subsample=FALSE, sequential=FALSE,verbose=FALSE,
                            mainClusterArgs=list(clusterArgs=list(evalClusterMethod="average")),
                            isCount=FALSE))
})

test_that("`getClusterManyParams` works", {
	cc<-clusterMany(mat, ks=c(3,4),nFilterDims=c(10,15),nPCADim=c(3,4),
		reduceMethod=c("none","PCA","var"),clusterFunction="pam",
	  	subsample=FALSE, sequential=FALSE,run=TRUE,verbose=FALSE,
	    isCount=FALSE)
	cc<-combineMany(cc,proportion=1,whichClusters = "clusterMany")
	expect_silent(paramAll<-getClusterManyParams(cc))
	expect_equal(colnames(paramAll),c("clusteringIndex", "reduceMethod", "nReducedDims", "nFilterDims", "k"))
	expect_true(is.data.frame(paramAll))
	expect_equal(sort(as.character(unique(paramAll[,"reduceMethod"]))),sort(c("none","var","PCA")))
	
	expect_equal(sort(unique(paramAll[,"nReducedDims"]),na.last=TRUE),sort(c(NA,3,4),na.last=TRUE))
	expect_equal(sort(unique(paramAll[,"nFilterDims"]),na.last=TRUE),sort(c(NA,10,15),na.last=TRUE))
	expect_true(is.numeric(paramAll[,"k"]))
	
	paramSub<-getClusterManyParams(cc,whichClusters=3:4)
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
    testSCF<-as(testSCE,"SingleCellFilter")
    
	#matrix
	expect_silent(ccVar<-clusterSingle(contData, 
	        subsample=FALSE, sequential=FALSE, reduceMethod="var",
	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	 	   isCount=FALSE))
   	expect_silent(ccPCA<-clusterSingle(contData, 
   	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
   	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
   	 	   isCount=FALSE))
  	expect_silent(ccNone<-clusterSingle(contData, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
  	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=FALSE))
	expect_silent(cm<-clusterMany(contData, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
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
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=FALSE))
  	expect_silent(ccPCA<-clusterSingle(testSE, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
  	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=FALSE))
 	expect_silent(ccNone<-clusterSingle(testSE, 
 	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
 	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	 	   isCount=FALSE))
	expect_silent(cm<-clusterMany(testSE, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
	   	        nReducedDims=3, nFilterDims=3,isCount=FALSE))
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 

	
   #SCE
	expect_silent(ccVar<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=FALSE))
  	expect_silent(ccPCA<-clusterSingle(testSCE, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
  	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=FALSE))
 	expect_silent(ccNone<-clusterSingle(testSCE, 
 	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
 	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	 	   isCount=FALSE))
	expect_silent(cm<-clusterMany(testSCE, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
	   	        nReducedDims=3, nFilterDims=3,isCount=FALSE))
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 
	
   #SCF
	expect_silent(ccVar<-clusterSingle(testSCF, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=FALSE))
  	expect_silent(ccPCA<-clusterSingle(testSCF, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
  	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=FALSE))
 	expect_silent(ccNone<-clusterSingle(testSCF, 
 	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
 	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	 	   isCount=FALSE))
	expect_silent(cm<-clusterMany(testSCF, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
	   	        nReducedDims=3, nFilterDims=3,isCount=FALSE))
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
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
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
    testSCF<-as(testSCE,"SingleCellFilter")
	expectTrans1<-round(log2(countData[1,]+1),2)

	#matrix
	expect_silent(ccVar<-clusterSingle(countData, 
	        subsample=FALSE, sequential=FALSE, reduceMethod="var",
	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
	 	   isCount=TRUE))
   	expect_silent(ccPCA<-clusterSingle(countData, 
   	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
   	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
   	 	   isCount=TRUE))
  	expect_silent(ccNone<-clusterSingle(countData, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
  	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=TRUE))
	expect_silent(cm<-clusterMany(countData, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
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
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
  	expect_silent(ccPCA<-clusterSingle(testSE, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
  	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=TRUE))
 	expect_silent(ccNone<-clusterSingle(testSE, 
 	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
 	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	 	   isCount=TRUE))
	expect_silent(cm<-clusterMany(testSE, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
	   	        nReducedDims=3, nFilterDims=3,isCount=TRUE))
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 

	
   #SCE
	expect_silent(ccVar<-clusterSingle(testSCE, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
  	expect_silent(ccPCA<-clusterSingle(testSCE, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
  	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=TRUE))
 	expect_silent(ccNone<-clusterSingle(testSCE, 
 	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
 	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	 	   isCount=TRUE))
	expect_silent(cm<-clusterMany(testSCE, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
	   	        nReducedDims=3, nFilterDims=3,isCount=TRUE))
	expect_equal(nClusterings(cm),3)	
	expect_silent(params<-getClusterManyParams(cm))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm)[1,],2), expectTrans1) 
	
   #SCF
	expect_silent(ccVar<-clusterSingle(testSCF, 
        subsample=FALSE, sequential=FALSE, reduceMethod="var",
        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	   isCount=TRUE))
  	expect_silent(ccPCA<-clusterSingle(testSCF, 
  	        subsample=FALSE, sequential=FALSE, reduceMethod="PCA",
  	        nDims=3, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
  	 	   isCount=TRUE))
 	expect_silent(ccNone<-clusterSingle(testSCF, 
 	        subsample=FALSE, sequential=FALSE, reduceMethod="none",
 	        nDims=NA, mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3)),
 	 	   isCount=TRUE))
	expect_silent(cm<-clusterMany(testSCF, clusterFunction="pam",ks=3,
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
	   	        nReducedDims=3, nFilterDims=3,isCount=TRUE))
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
	   	        subsample=FALSE, sequential=FALSE, reduceMethod=c("PCA","var","none"),
	   	        nReducedDims=3, nFilterDims=3))
	expect_equal(nClusterings(cm2),6)	
	expect_silent(params<-getClusterManyParams(cm2))	
   
	expect_equal(primaryCluster(ccPCA),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="PCA"]])
	expect_equal(primaryCluster(ccVar),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="var"]])
	expect_equal(primaryCluster(ccNone),clusterMatrix(cm2)[,params$clusteringIndex[params$reduceMethod=="none"]])
	expect_equal(round(transformData(cm2)[1,],2), expectTrans1) 
    
  
})
