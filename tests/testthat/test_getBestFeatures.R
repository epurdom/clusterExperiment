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

  ## add some unclustered
  top1 <- getBestFeatures(simData, primaryCluster(ceSimCont), contrastType="F",
                        isCount=FALSE)
  idx <- top1$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimCont)>0]), top1$AveExpr)

  ## check defaults
  topC0 <- getBestFeatures(ceSimCont)
  topC1 <- getBestFeatures(ceSimCont, contrastType="F",  isCount=FALSE)
  expect_equal(topC1, topC0)

  expect_equal(topC1, top1)

  top2 <- getBestFeatures(simData, primaryCluster(ceSimCont), contrastType="Pairs",
                        isCount=FALSE)
  idx <- top2$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimCont)>0]), top2$AveExpr)
  topC2 <- getBestFeatures(ceSimCont, contrastType="Pairs", isCount=FALSE)
  expect_equal(topC2, top2)

  top3 <- getBestFeatures(simData, primaryCluster(ceSimCont), contrastType="OneAgainstAll",
                        isCount=FALSE)
  idx <- top3$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimCont)>0]), top3$AveExpr)
  topC3 <- getBestFeatures(ceSimCont, contrastType="OneAgainstAll", 
                        isCount=FALSE)
  expect_equal(topC3, top3)

  ### test voom

  logcpm <- t(log2(t(simCount + 0.5)/(colSums(simCount) + 1) * 1e+06))
  voom1 <- getBestFeatures(simCount, primaryCluster(ceSim), contrastType="F",
                        isCount=TRUE)
  idx <- voom1$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom1$AveExpr)

  voom2 <- getBestFeatures(simCount, primaryCluster(ceSim), contrastType="Pairs",
                        isCount=TRUE)
  idx <- voom2$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom2$AveExpr)

  voom3 <- getBestFeatures(simCount, primaryCluster(ceSim), contrastType="OneAgainstAll",
                        isCount=TRUE)
  idx <- voom3$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom3$AveExpr)

  ## test dendrogram
  expect_error(getBestFeatures(simData, primaryCluster(ceSim), contrastType="Dendro"),
               "must provide dendro")
  
  dendro <- makeDendrogram(simData, primaryCluster(ceSimCont))
  expect_error(getBestFeatures(simData, primaryCluster(ceSimCont), contrastType="Dendro",
                               dendro=dendro$samples), "dendro don't match")
  

}
)
test_that("'Dendro' contrasts works for clusterExperiment object in `getBestFeatures`",{
  dendro <- makeDendrogram(simData, primaryCluster(ceSimCont))
  dend1 <- getBestFeatures(simData, primaryCluster(ceSimCont), contrastType="Dendro",
                           dendro = dendro$clusters)
  ceTemp<-ceSimCont
  ceTemp <- makeDendrogram(ceTemp)
  dendC1 <- getBestFeatures(ceTemp, contrastType="Dendro")
  expect_equal(dend1, dendC1)
  
  #check whole mergeDendrogram thing at least runs!
  cl1 <- clusterSingle(smSimData, 
                       subsample=FALSE, sequential=FALSE,
                       mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),isCount=FALSE)
  test<-clusterLegend(cl1)[[1]]
  test[,"name"]<-test[,"clusterIds"]
  clusterLegend(cl1)[[1]]<-test
  clustWithDendro <- makeDendrogram(cl1)
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP",plotInfo="none",plot=FALSE,calculateAll=FALSE)
  resCE<-getBestFeatures(clustMerged, contrastType="Dendro")

})

test_that("`plotContrastHeatmap` works", {
    topC2 <- getBestFeatures(ceSimCont, contrastType="Pairs", isCount=FALSE)
	plotContrastHeatmap(ceSimCont,signifTable=topC2)

    dendro <- makeDendrogram(ceSimCont, whichCluster=primaryClusterIndex(ceSimCont))
    topCD <- getBestFeatures(dendro, contrastType="Dendro", isCount=FALSE)
	plotContrastHeatmap(dendro,signifTable=topCD)
	
    top1 <- getBestFeatures(simData, primaryCluster(ceSimCont), contrastType="F",
                          isCount=FALSE)
	expect_error(plotContrastHeatmap(dendro,signifTable=top1),"signifTable must have columns 'IndexInOriginal' and 'Contrast'")
						  
})



