mat <- matrix(data=rnorm(200), ncol=10)
labels <- gl(5, 2)
se <- SummarizedExperiment(se, colData=DataFrame(clusterLabels=labels))

cc <- clusterCells(mat, labels, isLog = TRUE)
cc2 <- clusterCells(se, labels, isLog = TRUE)

stopifnot(identical(assay(clusterCells(assay(se), labels, isLog = TRUE)), assay(clusterCells(se, labels, isLog = TRUE))))
