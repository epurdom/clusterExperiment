context("InternalFunctions")
source("create_objects.R")

test_that("`.TypeIntoIndices` works",{
    ccM<-combineMany(cc,proportion=.2,clusterLabel="mySpecialLabel")
    ccM<-combineMany(ccM,proportion=.2)
    indCM<-clusterExperiment:::.TypeIntoIndices(ccM,wh=c("clusterMany"))
    expect_equal(length(indCM),nClusters(cc))
    expect_equal(clusterExperiment:::.TypeIntoIndices(ccM,wh=c("mySpecialLabel","combineMany")),c(2,1))
    expect_equal(clusterExperiment:::.TypeIntoIndices(ccM,wh=c("mySpecialLabel","clusterMany")),c(2,indCM))
})

