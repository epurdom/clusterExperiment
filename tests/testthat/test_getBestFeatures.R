context("getBestFeatures")
source("create_objects.R")

library(limma)
test_that("`clusterContrasts` works with matrix and ClusterExperiment objects", {
   x1<- clusterContrasts(primaryCluster(ccSE),contrastType="Pairs")
  x2<-clusterContrasts(ccSE,contrastType="Pairs")
  expect_equal(x1,x2)
  x1<- clusterContrasts(primaryCluster(ccSE),contrastType="OneAgainstAll")
  x2<-clusterContrasts(ccSE,contrastType="OneAgainstAll")
  expect_equal(x1,x2)
  expect_error(clusterContrasts(primaryCluster(ccSE),contrastType="Dendro"),"must provide dendrogram if contrastType='Dendro'")
  ccSE<-makeDendrogram(ccSE)
  x1<- clusterContrasts(primaryCluster(ccSE),contrastType="Dendro",dendro=ccSE@dendro_clusters)
  x2<-clusterContrasts(ccSE,contrastType="Dendro")
  expect_equal(x1,x2)
})
test_that("`getBestFeatures` works with matrix and ClusterExperiment objects", {

  cl <- clusterMany(simData, ks=c(3,4),clusterFunction="pam",
                              subsample=FALSE, sequential=FALSE,
                              isCount=FALSE,verbose=FALSE)
  ## add some unclustered
  clusters <- primaryCluster(cl)
  clusters[1:5] <- -1
  cl <- addClusters(cl, clusters)
  primaryClusterIndex(cl) <- 3
  table(primaryCluster(cl))

  top1 <- getBestFeatures(simData, primaryCluster(cl), contrastType="F",
                       returnType="Table", isCount=FALSE)
  idx <- top1$IndexInOriginal
  expect_equal(rowMeans(simData[idx,clusters>0]), top1$AveExpr)

  ## check defaults
  topC0 <- getBestFeatures(cl)
  topC1 <- getBestFeatures(cl, contrastType="F", returnType="Table", isCount=FALSE)
  expect_equal(topC1, topC0)

  expect_equal(topC1, top1)

  top2 <- getBestFeatures(simData, primaryCluster(cl), contrastType="Pairs",
                       returnType="Table", isCount=FALSE)
  idx <- top2$IndexInOriginal
  expect_equal(rowMeans(simData[idx,clusters>0]), top2$AveExpr)

  topC2 <- getBestFeatures(cl, contrastType="Pairs", returnType="Table",
                        isCount=FALSE)
  expect_equal(topC2, top2)

  top3 <- getBestFeatures(simData, primaryCluster(cl), contrastType="OneAgainstAll",
                       returnType="Table", isCount=FALSE)
  idx <- top3$IndexInOriginal
  expect_equal(rowMeans(simData[idx,clusters>0]), top3$AveExpr)

  topC3 <- getBestFeatures(cl, contrastType="OneAgainstAll", returnType="Table",
                        isCount=FALSE)
  expect_equal(topC3, top3)

  ### test voom

  logcpm <- t(log2(t(simCount + 0.5)/(colSums(simCount) + 1) * 1e+06))
  voom1 <- getBestFeatures(simCount, primaryCluster(cl), contrastType="F",
                       returnType="Table", isCount=TRUE)
  idx <- voom1$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,clusters>0]), voom1$AveExpr)

  voom2 <- getBestFeatures(simCount, primaryCluster(cl), contrastType="Pairs",
                       returnType="Table", isCount=TRUE)
  idx <- voom2$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,clusters>0]), voom2$AveExpr)

  voom3 <- getBestFeatures(simCount, primaryCluster(cl), contrastType="OneAgainstAll",
                       returnType="Table", isCount=TRUE)
  idx <- voom3$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,clusters>0]), voom3$AveExpr)

  ## test index output
  idx1 <- getBestFeatures(simData, primaryCluster(cl), contrastType="F",
                       returnType="Index", isCount=FALSE)
  expect_equal(as.numeric(idx1), top1$IndexInOriginal)

  idx2 <- getBestFeatures(simData, primaryCluster(cl), contrastType="Pairs",
                       returnType="Index", isCount=FALSE)
  expect_equal(as.numeric(idx2), top2$IndexInOriginal)

  idx3 <- getBestFeatures(simData, primaryCluster(cl), contrastType="OneAgainstAll",
                       returnType="Index", isCount=FALSE)
  expect_equal(as.numeric(idx3), top3$IndexInOriginal)

  ## test dendrogram
  expect_error(getBestFeatures(simData, primaryCluster(cl), contrastType="Dendro"),
               "must provide dendro")

  dendro <- makeDendrogram(simData, primaryCluster(cl))
  expect_error(getBestFeatures(simData, primaryCluster(cl), contrastType="Dendro",
                            dendro=dendro$samples), "dendro don't match")

  dend1 <- getBestFeatures(simData, primaryCluster(cl), contrastType="Dendro",
                        dendro = dendro$clusters)
  cl <- makeDendrogram(cl)
  dendC1 <- getBestFeatures(cl, contrastType="Dendro")
  expect_equal(dend1, dendC1)
}
)
