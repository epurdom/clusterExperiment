context("getBestFeatures")

plotAll<-FALSE #set to true to actually SEE the plots; otherwise for TravisCI, where no visual, runs quicker with FALSE
###Note some are still run with plot=TRUE to check that works with aheatmap. Only a fraction not plotted.
test_that("`clusterContrasts` works with matrix and ClusterExperiment objects", {
  expect_silent(x1<- clusterContrasts(primaryCluster(ccSE),contrastType="Pairs"))
  expect_silent(x2<-clusterContrasts(ccSE,contrastType="Pairs"))
  expect_equal(x1,x2)
  expect_silent(x1<- clusterContrasts(primaryCluster(ccSE),contrastType="OneAgainstAll"))
  expect_silent(x2<-clusterContrasts(ccSE,contrastType="OneAgainstAll"))
  expect_equal(x1,x2)
  expect_error(clusterContrasts(primaryCluster(ccSE),contrastType="Dendro"),"must provide dendrogram if contrastType='Dendro'")
  expect_silent(ccSE<-makeDendrogram(ccSE))
  expect_silent(x1<- clusterContrasts(primaryCluster(ccSE),contrastType="Dendro",dendro=ccSE@dendro_clusters))
  expect_false(any(is.na(x1$contrastNames)))
  expect_silent(x2<-clusterContrasts(ccSE,contrastType="Dendro"))
  expect_equal(x1,x2)
})



test_that("`getBestFeatures` works with matrix objects", {

  ## add some unclustered
  expect_silent(top1 <- getBestFeatures(simData, 
	  primaryCluster(ceSimData), contrastType="F", DEMethod="limma")
						)
  idx <- top1$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimData)>0]), top1$AveExpr)


  expect_silent(top2 <- getBestFeatures(simData, 
	  primaryCluster(ceSimData), contrastType="Pairs",
                        DEMethod="limma"))
  idx <- top2$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimData)>0]), top2$AveExpr)
  expect_silent(topC2 <- getBestFeatures(ceSimData, contrastType="Pairs",DEMethod="limma"))
  expect_equal(topC2, top2)

  expect_silent(top3 <- getBestFeatures(simData, 
	  primaryCluster(ceSimData), contrastType="OneAgainstAll",
                        DEMethod="limma"))
  idx <- top3$IndexInOriginal
  expect_equal(rowMeans(simData[idx,primaryCluster(ceSimData)>0]), top3$AveExpr)
  expect_silent(topC3 <- getBestFeatures(ceSimData, 
	  contrastType="OneAgainstAll", 
      DEMethod="limma"))
  expect_equal(topC3, top3)

  ### test voom

  logcpm <- t(log2(t(simCount + 0.5)/(colSums(simCount) + 1) * 1e+06))
  expect_silent(voom1 <- getBestFeatures(simCount, 
	  primaryCluster(ceSim), contrastType="F",
                        DEMethod="limma-voom"))
  idx <- voom1$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom1$AveExpr)

  expect_silent(voom2 <- getBestFeatures(simCount,
	   primaryCluster(ceSim), contrastType="Pairs",
                        DEMethod="limma-voom"))
  idx <- voom2$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom2$AveExpr)

  expect_silent(voom3 <- getBestFeatures(simCount, primaryCluster(ceSim), contrastType="OneAgainstAll",
                        DEMethod="limma-voom"))
  idx <- voom3$IndexInOriginal
  expect_equal(rowMeans(logcpm[idx,primaryCluster(ceSim)>0]), voom3$AveExpr)

}
)

test_that("`getBestFeatures` works with ClusterExperiment objects", {
    ## check F
    expect_silent(top1 <- getBestFeatures(simData, 
  	  primaryCluster(ceSimData), contrastType="F", DEMethod="limma")
    expect_silent(topC0 <- getBestFeatures(ceSimData, DEMethod="limma"))
    expect_silent(topC1 <- getBestFeatures(ceSimData, contrastType="F",  DEMethod="limma"))
    expect_equal(topC1, topC0)
    expect_equal(topC1, top1)
	
	##check return format of contrasts
	# > clusterLegend(ceSimData)[[primaryClusterIndex(ceSimData)]]
# 	     clusterIds color     name
# 	[1,] "-1"       "white"   "-1"
# 	[2,] "-2"       "grey"    "-2"
# 	[3,] "1"        "#E31A1C" "a"
# 	[4,] "2"        "#1F78B4" "b"
# 	[5,] "3"        "#33A02C" "c"
# 	[6,] "4"        "#FF7F00" "d"
    expect_silent(topPairs <- getBestFeatures(ceSimData, DEMethod="limma",contrastType="Pairs"))
    expect_true(topPairs$ContrastName[1]=="a-b")
    expect_true(topPairs$InternalName[1]=="Cl01-Cl02")

    expect_silent(topOne <- getBestFeatures(ceSimData, DEMethod="limma",contrastType="OneAgainstAll"))
    expect_true(topOne$ContrastName[1]=="a")
    expect_true(topOne$InternalName[1]=="Cl01")
	
	expect_silent(ceDend <- makeDendrogram(ceSimData))
	expect_silent(topDend<-getBestFeatures(ceDend, DEMethod="limma",contrastType="Dendro"))
	##Need to add correct expection after work out labels.
})
test_that("'Dendro' contrasts works for ClusterExperiment object in `getBestFeatures`",{
  ## test dendrogram
  expect_error(getBestFeatures(simData, primaryCluster(ceSim), contrastType="Dendro"),
               "must provide dendro")
  
  expect_silent(dendro <- makeDendrogram(simData, primaryCluster(ceSimData)))
  expect_error(getBestFeatures(simData, primaryCluster(ceSimData), contrastType="Dendro",
                               dendro=dendro$samples,DEMethod="limma"), "dendro don't match")
  expect_silent(dendro <- makeDendrogram(simData, primaryCluster(ceSimData)))
  expect_silent(dend1 <- getBestFeatures(simData, primaryCluster(ceSimData), contrastType="Dendro", dendro = dendro$clusters,DEMethod="limma"))
						   
  length(grep("NodeId",dend1$ContrastName))
  ceTemp<-ceSimData
  expect_silent(ceTemp <- makeDendrogram(ceTemp))
  expect_silent(dendC1 <- getBestFeatures(ceTemp, contrastType="Dendro",DEMethod="limma"))
  expect_equal(dend1, dendC1)
  
  #check whole mergeDendrogram thing at least runs!
  expect_silent(cl1 <- clusterSingle(smSimData, 
                       subsample=FALSE, sequential=FALSE,
                       mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)), 
											 isCount=FALSE)
											 )
  expect_silent(clustWithDendro <- makeDendrogram(cl1))
  expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP",plotInfo="none",plot=FALSE,calculateAll=FALSE, DEMethod="limma"),"Merging will be done on")
  expect_silent(getBestFeatures(clustMerged, contrastType="Dendro",DEMethod="limma"))

})
test_that("`getBestFeatures` works with HDF5 assay slot",{
    expect_silent(cl1 <- clusterSingle(hdfObj, 
            subsample=FALSE, sequential=FALSE,
			mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
			isCount=FALSE))
    expect_silent(getBestFeatures(cl1,DEMethod="limma"))
								
	
})
test_that("`plotContrastHeatmap` works", {
	expect_silent(ceSimData<-renameClusters(ceSimData,whichCluster=1,val=letters[1:nClusters(ceSimData)[1]]))
    expect_silent(topC2 <- getBestFeatures(ceSimData, contrastType="Pairs", DEMethod="limma"))
	expect_silent(plotContrastHeatmap(ceSimData,signifTable=topC2))

	expect_silent(topCOne <- getBestFeatures(ceSimData, contrastType="OneAgainstAll", DEMethod="limma"))
	expect_silent(plotContrastHeatmap(ceSimData,signifTable=topCOne,plot=plotAll))
	  
    dendro <- makeDendrogram(ceSimData, whichCluster=primaryClusterIndex(ceSimData))
	# > 	phylobase::nodeLabels(dendro@dendro_clusters)
# 	      5       6       7
# 	"Node1" "Node2" "Node3"
	topCD <- getBestFeatures(dendro, contrastType="Dendro", DEMethod="limma")
	plotContrastHeatmap(dendro,signifTable=topCD,plot=plotAll)
	
    top1 <- getBestFeatures(simData, primaryCluster(ceSimData), contrastType="F",
                          DEMethod="limma")
	expect_error(plotContrastHeatmap(dendro,signifTable=top1),"signifTable must have columns 'IndexInOriginal' and 'Contrast'")
						  
	#test name replacement:
	plotContrastHeatmap(ceSimData,signifTable=topC2,whichCluster=primaryClusterIndex(ceSimData),plot=plotAll)
	plotContrastHeatmap(ceSimData,signifTable=topCOne,whichCluster=primaryClusterIndex(ceSimData),plot=plotAll)
	plotContrastHeatmap(ceSimData,signifTable=topCD,whichCluster=primaryClusterIndex(ceSimData),plot=plotAll)
	expect_error(plotContrastHeatmap(ceSimData,signifTable=topC2,whichCluster=c(1,2)),"must identify only a single clustering")
	expect_error(plotContrastHeatmap(ceSimData,signifTable=topC2,whichCluster=50),"Invalid value for 'whichCluster'. Must be integer between")
	
})

test_that("`getBestFeatures` works with weights", {
	set.seed(258179)
	weights<-matrix(runif(nrow(seSimCount)*ncol(seSimCount)),nrow=nrow(seSimCount),ncol=ncol(seSimCount))
	ceSimW<-ceSim
	assay(ceSimW,"weights")<-weights

	expect_silent(outW<-getBestFeatures(ceSim,weights=weights,DEMethod="edgeR",contrastType="F"))
	expect_silent(outNW<-getBestFeatures(ceSim,weights=NULL,DEMethod="edgeR",contrastType="F"))
	expect_silent(outV<-getBestFeatures(ceSim,weights=weights,DEMethod="limma-voom",contrastType="F"))
	expect_silent(outV2<-getBestFeatures(ceSim,DEMethod="limma-voom",contrastType="F"))
	expect_equal(outV,outV2)

	expect_silent(outNW_C<-getBestFeatures(ceSim,weights=NULL,DEMethod="edgeR",contrastType="Pairs"))
	expect_silent(outW_C<-getBestFeatures(ceSim,weights=weights,DEMethod="edgeR",contrastType="Pairs"))
	expect_silent(getBestFeatures(ceSim,weights=weights,DEMethod="edgeR",contrastType="Pairs",contrastAdj="PerContrast"))
	expect_silent(getBestFeatures(ceSim,weights=weights,DEMethod="edgeR",contrastType="Pairs",contrastAdj="AfterF"))
	
	#-------
	##Test get right answers with weights (have to remove unassigned)
	#-------
	expect_silent(cleanSim<-removeUnassigned(ceSimW,whichCluster="primary"))
	expect_true(all(primaryCluster(cleanSim)>0))
	expect_silent(pairsContrasts<-clusterContrasts(cleanSim,contrastType="Pairs"))
  designContr <- model.matrix(~ 0 + factor(primaryCluster(cleanSim)))
  fitContr<- edgeR::DGEList(assay(cleanSim,1))
  fitContr$weights <- assay(cleanSim,"weights")
  fitContr <- edgeR::estimateDisp(fitContr, designContr)
  fitContr <- edgeR::glmFit(fitContr,design = designContr)
  lrt <- zinbwave::glmWeightedF(fitContr, contrast = pairsContrasts$contrastMatrix[,1,drop=FALSE])
  resultsManual<- edgeR::topTags(lrt)
	expect_silent(outWClean<-getBestFeatures(cleanSim,DEMethod="edgeR",contrastAdj="PerContrast",contrastType="Pairs"))
	resultsC1<-outWClean[outWClean$Contrast == colnames(pairsContrasts$contrastMatrix)[1],]
	expect_equal(unname(data.matrix(resultsC1[,c("logFC","logCPM","LR","P.Value","padjFilter","adj.P.Val")])), unname(data.matrix(data.frame(resultsManual))))
	
	### check when weights in assay picks them up automatically (when have right name...)
	expect_silent(outNW_C2<-getBestFeatures(ceSimW,weights=NULL,DEMethod="edgeR",contrastType="Pairs"))
	expect_equal(outNW_C2,outNW_C)

	expect_silent(outW_C2<-getBestFeatures(ceSimW,DEMethod="edgeR",contrastType="Pairs"))
	expect_equal(outW_C2,outW_C)
	
	##Test mergeClusters
	expect_silent(ceSimW<-makeDendrogram(ceSimW))
	expect_message(mergeClusters(ceSimW,DEMethod="edgeR"),"Merging will be done on")
	expect_message(mergeClusters(ceSimW,DEMethod="edgeR",weights=NULL,forceCalculate=TRUE),"Merging will be done on")
	
	##Test RSEC
	expect_message(RSEC(ceSimW,sequential=FALSE,subsample=FALSE,ks=4:15,isCount=TRUE,stopOnError=TRUE),"merging with these parameters did not result in any clusters being merged")
	
	sceCountData<-SingleCellExperiment(simCount)
	reducedDims(sceCountData) <- reducedDims(sceSimDataDimRed)	
	expect_message(RSEC(sceCountData,sequential=TRUE,subsample=FALSE,clusterFunction="pam",k0s=4:15,isCount=TRUE, stopOnError=TRUE),"Merging will be done on")
	
})
