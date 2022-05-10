test_that("getReducedData works as with clusterSingle",{
  expect_silent(clustNothingDimRed <- 
    clusterSingle(sceSimDataDimRed,
       reduceMethod = "none",
       mainClusterArgs=list( clusterArgs=list(k=3),
       clusterFunction=listBuiltInFunctions()[[1]]),
       subsample=FALSE, sequential=FALSE, isCount=FALSE))  	
  expect_silent(clustNothing <- 
    clusterSingle(sceSimData,
      reduceMethod = "none",
      mainClusterArgs=list( clusterArgs=list(k=3),
      clusterFunction=listBuiltInFunctions()[[1]]),
      subsample=FALSE, sequential=FALSE, isCount=FALSE))  	

  #dimReduce
  expect_silent(out1<-getReducedData(clustNothing,reduceMethod="PCA"))
  expect_true("PCA" %in% reducedDimNames(out1))
  expect_equal(colData(out1),colData(clustNothing))
  expect_equal(rownames(out1),rownames(clustNothing))
  expect_equal(colnames(out1),colnames(clustNothing))
  expect_equal(metadata(out1),metadata(clustNothing))
  
  
  expect_warning(getReducedData(clustNothingDimRed,reduceMethod="PCA"),
    "will not add reduced dataset to object because already exists method with that name")
  ## FIXME:
  #expect_silent(
  out1<-getReducedData(clustNothing,
      reduceMethod="PCA",reducedDimName="MyPCA")
  #)
  expect_true("MyPCA" %in% reducedDimNames(out1))
  expect_false("PCA" %in% reducedDimNames(out1))
  
  expect_warning(getReducedData(clustNothingDimRed,
      reduceMethod="PCA",nDim=200,reducedDimName="MyPCA",
      returnValue="list"),
      "requesting an existing dimensionality reduction but with greater number of dimensions than available")
      expect_silent(out3<-getReducedData(clustNothingDimRed,reduceMethod="PCA",nDim=50,reducedDimName="MyPCA",returnValue="list"))
  expect_equal(nrow(out3$dataMatrix),50)
  
  #filters
  expect_silent(out2<-getReducedData(clustNothing,reduceMethod="mad",nDim=50))
  expect_true("mad_clusterSingle" %in% filterNames(out2))
  expect_true("filteredBy_mad_clusterSingle" %in% reducedDimNames(out2))
  expect_equal(ncol(reducedDim(out2,"filteredBy_mad_clusterSingle")),50)
  
  expect_silent(out5<-getReducedData(clustNothing,reduceMethod="mad",filterIgnoresUnassigned=FALSE))
  expect_true("mad" %in% filterNames(out5))
  expect_true("filteredBy_mad" %in% reducedDimNames(out5))
  
})
