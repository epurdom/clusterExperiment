context("combineMany")
source("create_objects.R")

test_that("`combineMany` works with matrix and ClusterExperiment objects", {
            clustNothing <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
                                        subsample=FALSE, sequential=FALSE,
                                        isCount=FALSE,verbose=FALSE)
            x1<-combineMany(clustNothing,proportion=1,whichClusters = "clusterMany")
            x2<-combineMany(clustNothing,proportion=1)
            expect_equal(x1,x2)
			ceObj<-clusterSingle(mat, subsample=FALSE,
				clusterDArgs=list(clusterFunction="pam",clusterArgs=list(k=3)))
            expect_error(combineMany(ceObj,proportion=1),
                         "no clusters specified")
			expect_error(combineMany(ceObj,whichCluster="clusterSingle"),
			                          'argument "proportion" is missing, with no default')

            shared1 <- combineMany(clusterMatrix(clustNothing),proportion=1)
            shared2 <- combineMany(clustNothing, "all",proportion=1)
            expect_equal(shared1$clustering, primaryCluster(shared2))

            shared3 <- combineMany(clustNothing, "workflow",proportion=1)
            expect_equal(shared2, shared3)

            shared4 <- combineMany(clustNothing, 1:nClusters(clustNothing),proportion=1)
            expect_equal(shared3, shared4)

            shared5 <- combineMany(clustNothing, "workflow",
                                          proportion=.5)

            shared6 <- combineMany(clustNothing, "workflow",
                                          proportion=.5,
                                          clusterFunction="tight")

            expect_error(combineMany(clustNothing, "workflow",
                                          proportion=.5,
                                          clusterFunction="pam"), "implemented")

            expect_true("combineMany" %in% clusterTypes(shared6))
            expect_true("combineMany" %in% clusterLabels(shared6))
            expect_true(all(primaryCluster(shared6)==clusterMatrix(shared6)[,"combineMany"]))
})

test_that("`combineMany` works when multiple runs of workflow", {
  clustNothing <- clusterMany(mat, ks=c(3,4),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE)

  shared1 <- combineMany(clustNothing, "all",proportion=1)
  shared2<-combineMany(shared1,"all",proportion=1)
  expect_true("combineMany.1" %in% clusterTypes(shared2))

  clustNothing2 <- clusterMany(shared2, ks=c(5,6), clusterFunction="pam",
                               subsample=FALSE, sequential=FALSE,
                               isCount=FALSE,verbose=FALSE)
  expect_true("combineMany.1" %in% clusterTypes(clustNothing2))
  expect_true("clusterMany.2" %in% clusterTypes(clustNothing2))
  expect_true("combineMany.2" %in% clusterTypes(clustNothing2))

  shared3 <- combineMany(clustNothing2, "all",proportion=1)
  shared4 <- combineMany(clusterMatrix(clustNothing2),proportion=1)
  expect_equal(shared4$clustering, primaryCluster(shared3))

  shared5 <- combineMany(clustNothing2, "workflow",proportion=1)
  shared6 <- combineMany(clusterMatrix(clustNothing2)[,1:2],proportion=1)
  expect_equal(shared6$clustering, primaryCluster(shared5))


  clustNothing3 <- addClusters(clustNothing2, primaryCluster(shared5))
  shared7 <- combineMany(clustNothing3, "all",proportion=1)
  shared8 <- combineMany(clustNothing3, "workflow",proportion=1)
})

test_that("`combineMany` preserves the colData and rowData of SE", {
  cl <- clusterMany(se, clusterFunction="pam", ks=2:4, findBestK=c(FALSE),
                    removeSil=TRUE, subsample=FALSE)

  cl <- combineMany(cl, whichClusters="workflow", proportion=0.5)

  expect_equal(colData(cl),colData(se))
  expect_equal(rownames(cl),rownames(se))
  expect_equal(colnames(cl),colnames(se))
  expect_equal(metadata(cl),metadata(se))
  expect_equal(rowData(cl),rowData(se))

})
