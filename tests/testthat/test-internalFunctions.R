context("InternalFunctions")


test_that("`.TypeIntoIndices` works",{
    ccM<-makeConsensus(ceSim,proportion=.2,clusterLabel="mySpecialLabel")
    ccM<-makeConsensus(ccM,proportion=.2)
    indCM<-clusterExperiment:::.TypeIntoIndices(ccM,wh=c("clusterMany"))
    expect_equal(length(indCM),sum(clusterTypes(ccM)=="clusterMany"))
    expect_equal(clusterExperiment:::.TypeIntoIndices(ccM,wh=c("mySpecialLabel","makeConsensus")),c(2,1))
    expect_equal(clusterExperiment:::.TypeIntoIndices(ccM,wh=c("mySpecialLabel","clusterMany")),c(2,indCM))
})

