library(clusterExperiment)
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version. Can delete this once package is stabilized.

test_that("`getBestGenes` works with matrix and ClusterExperiment objects", {

  cl <- clusterMany(simData, ks=c(3,4),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE)
  ## add some unclustered
  clusters <- primaryCluster(cl)
  clusters[1:5] <- -1
  cl <- addClusters(cl, clusters)
  primaryClusterIndex(cl) <- 3
  table(primaryCluster(cl))

  top1 <- getBestGenes(simData, primaryCluster(cl), type="F",
                       returnType="Table", voomCorrection=FALSE)
  idx <- top1$IndexInOriginal
  expect_equal(rowMeans(simData[idx,clusters>0]), top1$AveExpr)

  ## check defaults
  topC0 <- getBestGenes(cl)
  topC1 <- getBestGenes(cl, type="F", returnType="Table", voomCorrection=FALSE)
  expect_equal(topC1, topC0)

  expect_equal(topC1, top1)

  top2 <- getBestGenes(simData, primaryCluster(cl), type="Pairs",
                       returnType="Table", voomCorrection=FALSE)
  idx <- top2$IndexInOriginal
  expect_equal(rowMeans(simData[idx,clusters>0]), top2$AveExpr)

  topC2 <- getBestGenes(cl, type="Pairs", returnType="Table",
                        voomCorrection=FALSE)
  expect_equal(topC2, top2)

  top3 <- getBestGenes(simData, primaryCluster(cl), type="OneAgainstAll",
                       returnType="Table", voomCorrection=FALSE)
  idx <- top3$IndexInOriginal
  expect_equal(rowMeans(simData[idx,clusters>0]), top3$AveExpr)

  topC3 <- getBestGenes(cl, type="OneAgainstAll", returnType="Table",
                        voomCorrection=FALSE)
  expect_equal(topC3, top3)

  ### test voom

  logcpm <- t(log2(t(simCount + 0.5)/(colSums(simCount) + 1) * 1e+06))
  voom1 <- getBestGenes(simCount, primaryCluster(cl), type="F",
                       returnType="Table", voomCorrection=TRUE)
  idx <- voom1$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,clusters>0]), voom1$AveExpr)

  voom2 <- getBestGenes(simCount, primaryCluster(cl), type="Pairs",
                       returnType="Table", voomCorrection=TRUE)
  idx <- voom2$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,clusters>0]), voom2$AveExpr)

  voom3 <- getBestGenes(simCount, primaryCluster(cl), type="OneAgainstAll",
                       returnType="Table", voomCorrection=TRUE)
  idx <- voom3$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,clusters>0]), voom3$AveExpr)

  ## test index output
  idx1 <- getBestGenes(simData, primaryCluster(cl), type="F",
                       returnType="Index", voomCorrection=FALSE)
  expect_equal(as.numeric(idx1), top1$IndexInOriginal)

  idx2 <- getBestGenes(simData, primaryCluster(cl), type="Pairs",
                       returnType="Index", voomCorrection=FALSE)
  expect_equal(as.numeric(idx2), top2$IndexInOriginal)

  idx3 <- getBestGenes(simData, primaryCluster(cl), type="OneAgainstAll",
                       returnType="Index", voomCorrection=FALSE)
  expect_equal(as.numeric(idx3), top3$IndexInOriginal)

  ## test dendrogram
  expect_error(getBestGenes(simData, primaryCluster(cl), type="Dendro"),
               "must provide dendro")

  expect_error(getBestGenes(simData, primaryCluster(cl), type="Dendro",
                            dendro=dendrogram(cl)), "dendro don't match")

  dend1 <- getBestGenes(simData, primaryCluster(cl), type="Dendro",
                        dendro = makeDendrogram(cl, leaves="clusters"))
  dendC1 <- getBestGenes(cl, type="Dendro")
  expect_equal(dend1, dendC1)
}
)
