context("getBestFeatures")

plotAll<-FALSE #set to true to actually SEE the plots; otherwise for TravisCI, where no visual, runs quicker with FALSE
###Note some are still run with plot=TRUE to check that works with aheatmap. Only a fraction not plotted.
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


test_that("`getBestFeatures` works with HDF5 assay slot",{
    expect_silent(cl1 <- clusterSingle(hdfObj, 
            subsample=FALSE, sequential=FALSE,
			mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
			isCount=FALSE))
    expect_silent(getBestFeatures(cl1,isCount=FALSE))
								
	
})
test_that("`getBestFeatures` works with matrix and ClusterExperiment objects", {

  ## add some unclustered
  expect_silent(top1 <- getBestFeatures(simData, 
	  primaryCluster(ceSimData), contrastType="F",
                        isCount=FALSE)
						)
  idx <- top1$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimData)>0]), top1$AveExpr)

  ## check defaults
  expect_silent(topC0 <- getBestFeatures(ceSimData))
  expect_silent(topC1 <- getBestFeatures(ceSimData, contrastType="F",  isCount=FALSE))
  expect_equal(topC1, topC0)

  expect_equal(topC1, top1)

  expect_silent(top2 <- getBestFeatures(simData, 
	  primaryCluster(ceSimData), contrastType="Pairs",
                        isCount=FALSE))
  idx <- top2$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimData)>0]), top2$AveExpr)
  expect_silent(topC2 <- getBestFeatures(ceSimData, contrastType="Pairs", isCount=FALSE))
  expect_equal(topC2, top2)

  expect_silent(top3 <- getBestFeatures(simData, 
	  primaryCluster(ceSimData), contrastType="OneAgainstAll",
                        isCount=FALSE))
  idx <- top3$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimData)>0]), top3$AveExpr)
  expect_silent(topC3 <- getBestFeatures(ceSimData, 
	  contrastType="OneAgainstAll", 
      isCount=FALSE))
  expect_equal(topC3, top3)

  ### test voom

  logcpm <- t(log2(t(simCount + 0.5)/(colSums(simCount) + 1) * 1e+06))
  expect_silent(voom1 <- getBestFeatures(simCount, 
	  primaryCluster(ceSim), contrastType="F",
                        isCount=TRUE))
  idx <- voom1$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom1$AveExpr)

  expect_silent(voom2 <- getBestFeatures(simCount,
	   primaryCluster(ceSim), contrastType="Pairs",
                        isCount=TRUE))
  idx <- voom2$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom2$AveExpr)

  expect_silent(voom3 <- getBestFeatures(simCount, primaryCluster(ceSim), contrastType="OneAgainstAll",
                        isCount=TRUE))
  idx <- voom3$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom3$AveExpr)

  

}
)
test_that("'Dendro' contrasts works for ClusterExperiment object in `getBestFeatures`",{
  ## test dendrogram
  expect_error(getBestFeatures(simData, primaryCluster(ceSim), contrastType="Dendro"),
               "must provide dendro")
  
  dendro <- makeDendrogram(simData, primaryCluster(ceSimData))
  expect_error(getBestFeatures(simData, primaryCluster(ceSimData), contrastType="Dendro",
                               dendro=dendro$samples), "dendro don't match")
  dendro <- makeDendrogram(simData, primaryCluster(ceSimData))
  dend1 <- getBestFeatures(simData, primaryCluster(ceSimData), contrastType="Dendro",
                           dendro = dendro$clusters)
  ceTemp<-ceSimData
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
    mat<-clusterLegend(ceSimData)[[1]]
    mat[,"name"]<-letters[1:nrow(mat)]
    clusterLegend(ceSimData)[[1]]<-mat
    topC2 <- getBestFeatures(ceSimData, contrastType="Pairs", isCount=FALSE)
	  plotContrastHeatmap(ceSimData,signifTable=topC2)

	  topCOne <- getBestFeatures(ceSimData, contrastType="OneAgainstAll", isCount=FALSE)
	  plotContrastHeatmap(ceSimData,signifTable=topCOne,plot=plotAll)
	  
    dendro <- makeDendrogram(ceSimData, whichCluster=primaryClusterIndex(ceSimData))
    topCD <- getBestFeatures(dendro, contrastType="Dendro", isCount=FALSE)
	plotContrastHeatmap(dendro,signifTable=topCD,plot=plotAll)
	
    top1 <- getBestFeatures(simData, primaryCluster(ceSimData), contrastType="F",
                          isCount=FALSE)
	expect_error(plotContrastHeatmap(dendro,signifTable=top1),"signifTable must have columns 'IndexInOriginal' and 'Contrast'")
						  
	#test name replacement:
	plotContrastHeatmap(ceSimData,signifTable=topC2,whichCluster=primaryClusterIndex(ceSimData),plot=plotAll)
	plotContrastHeatmap(ceSimData,signifTable=topCOne,whichCluster=primaryClusterIndex(ceSimData),plot=plotAll)
	plotContrastHeatmap(ceSimData,signifTable=topCD,whichCluster=primaryClusterIndex(ceSimData),plot=plotAll)
	expect_error(plotContrastHeatmap(ceSimData,signifTable=topC2,whichCluster=c(1,2)),"Must indicate single clustering in 'whichCluster'")
	expect_error(plotContrastHeatmap(ceSimData,signifTable=topC2,whichCluster=50),"Did not indicate valid cluster in whichCluster argument")
	
})



