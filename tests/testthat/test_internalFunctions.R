context("InternalFunctions")


test_that("`.TypeIntoIndices` works",{
    ccM<-combineMany(ceSim,proportion=.2,clusterLabel="mySpecialLabel")
    ccM<-combineMany(ccM,proportion=.2)
    indCM<-clusterExperiment:::.TypeIntoIndices(ccM,wh=c("clusterMany"))
    expect_equal(length(indCM),sum(clusterTypes(ccM)=="clusterMany"))
    expect_equal(clusterExperiment:::.TypeIntoIndices(ccM,wh=c("mySpecialLabel","combineMany")),c(2,1))
    expect_equal(clusterExperiment:::.TypeIntoIndices(ccM,wh=c("mySpecialLabel","clusterMany")),c(2,indCM))
})

