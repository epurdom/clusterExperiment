context("combineMany")


test_that("`combineMany` works with matrix and ClusterExperiment objects", {
    expect_silent(clustNothing <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
                                subsample=FALSE, sequential=FALSE,
                                isCount=FALSE,verbose=FALSE))
    expect_silent(x1<-combineMany(clustNothing,proportion=1,whichClusters = "clusterMany"))
    expect_message(x2<-combineMany(clustNothing,proportion=1),"no clusters specified to combine")
    expect_equal(x1,x2)
	expect_silent(ceObj<-clusterSingle(mat, subsample=FALSE,
		mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=3))))


    expect_silent(shared1 <- combineMany(clusterMatrix(clustNothing),proportion=1))
    expect_silent(shared2 <- combineMany(clustNothing, "all",proportion=1))
    expect_equal(shared1$clustering, primaryCluster(shared2))

    expect_silent(shared3 <- combineMany(clustNothing, "workflow",proportion=1))
    expect_equal(shared2, shared3)

    expect_silent(shared4 <- combineMany(clustNothing, 1:nClusterings(clustNothing),proportion=1))
    expect_equal(shared3, shared4)

    expect_silent(shared5 <- combineMany(clustNothing, "workflow", proportion=.5))

    expect_silent(shared6 <- combineMany(clustNothing, "workflow",
                                  proportion=.5, clusterFunction="tight"))
	expect_true("combineMany" %in% clusterTypes(shared6))
	expect_true("combineMany" %in% clusterLabels(shared6))
	expect_true(all(primaryCluster(shared6)==clusterMatrix(shared6)[,"combineMany"]))

  	#----
  	##Check error catching:
  	#----
  	expect_error(combineMany(ceObj,proportion=1),
                   "no clusters specified")
 	expect_error(combineMany(ceObj,whichClusters = "all",proportion=-1),
           "Invalid value for the 'proportion' parameter")
	expect_error(combineMany(ceObj,whichClusters = "all",proportion=0.7, minSize=-3),
	       "Invalid value for the 'minSize' parameter")
	expect_error(combineMany(ceObj,whichClusters = "all",proportion=0.7, propUnassigned=3),
	       "Invalid value for the 'propUnassigned' parameter")
	expect_error(combineMany(ceObj,whichCluster="clusterSingle"),
  	        'argument "proportion" is missing, with no default')
  	expect_error(combineMany(ceObj,whichCluster="clusterSingle"),
  			'argument "proportion" is missing, with no default')
    expect_error(combineMany(clustNothing, "workflow",
           proportion=.5,clusterFunction="pam"), 
		   "combineMany is only implemented for '01' type")

})

test_that("`combineMany` works with hdf5",{
    expect_silent(clustNothing <- clusterMany(hdfSCE,
		 ks=c(3,4),clusterFunction="pam",
		 subsample=FALSE, sequential=FALSE,
		 isCount=FALSE,verbose=FALSE))
    expect_silent(x1<-combineMany(clustNothing,proportion=1,whichClusters = "clusterMany"))
	
})

test_that("`combineMany` works when multiple runs of workflow", {
  expect_silent(clustNothing <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE))

  expect_silent(shared1 <- combineMany(clustNothing, "all",proportion=1))
  expect_silent(shared2<-combineMany(shared1,"all",proportion=1))
  expect_true("combineMany.1" %in% clusterTypes(shared2))

  expect_silent(clustNothing2 <- clusterMany(shared2, ks=c(5,6), clusterFunction="pam",
                               subsample=FALSE, sequential=FALSE, verbose=FALSE))
  expect_true("combineMany.1" %in% clusterTypes(clustNothing2))
  expect_true("clusterMany.2" %in% clusterTypes(clustNothing2))
  expect_true("combineMany.2" %in% clusterTypes(clustNothing2))

  expect_silent(shared3 <- combineMany(clustNothing2, "all",proportion=1))
  expect_silent(shared4 <- combineMany(clusterMatrix(clustNothing2),proportion=1))
  expect_equal(shared4$clustering, primaryCluster(shared3))

  expect_silent(shared5 <- combineMany(clustNothing2, "workflow",proportion=1))
  expect_silent(shared6 <- combineMany(clusterMatrix(clustNothing2)[,1:2],proportion=1))
  expect_equal(shared6$clustering, primaryCluster(shared5))


  expect_silent(clustNothing3 <- addClusterings(clustNothing2, primaryCluster(shared5)))
  expect_silent(shared7 <- combineMany(clustNothing3, "all",proportion=1))
  expect_silent(shared8 <- combineMany(clustNothing3, "workflow",proportion=1))
})

test_that("`combineMany` preserves the colData and rowData of SE", {
  expect_silent(cl <- clusterMany(se, clusterFunction="pam", ks=2:4, findBestK=c(FALSE),
                    removeSil=TRUE, subsample=FALSE))

  expect_silent(cl <- combineMany(cl, whichClusters="workflow", proportion=0.5))

  expect_equal(colData(cl),colData(se))
  expect_equal(rownames(cl),rownames(se))
  expect_equal(colnames(cl),colnames(se))
  expect_equal(metadata(cl),metadata(se))
  expect_equal(rowData(cl),rowData(se))

})
