library(testthat)
library(clusterCells)

test_check("clusterCells")
expect_error(compareChoices(simData, ks=5,
                             alphas=c(0.1), findBestK=c(TRUE,FALSE),
                             sequential=c(FALSE),
                             subsample=c(TRUE),
                             removeSil=c(TRUE,FALSE), clusterMethod=c("pam","tight","hierarchical"),
                             clusterDArgs = list(minSize = 5,kRange=2:15),ncores=1,random.seed=48120)
             ,"must provide k in subsampleArgs for those with findBestK=TRUE and sequential=FALSE")
