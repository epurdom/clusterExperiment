context("getBestFeatures")
#check whole mergeDendrogram thing at least runs!
test_that("'Dendro' contrasts works for ClusterExperiment object in `getBestFeatures`",{
  expect_silent(cl1 <- clusterSingle(smSimData, 
      subsample=FALSE, sequential=FALSE,
      mainClusterArgs=list(clusterFunction="pam",
          clusterArgs=list(k=6)), 
			isCount=FALSE)
	)
  expect_silent(clustWithDendro <- makeDendrogram(cl1))
  expect_message(clustMerged <- mergeClusters(clustWithDendro, 
    mergeMethod="adjP",plotInfo="none",plot=FALSE,
    calculateAll=FALSE, DEMethod="limma"),"Merging will be done on")
  expect_silent(getBestFeatures(clustMerged, 
    contrastType="Dendro",DEMethod="limma"))
  })
  
test_that("`getBestFeatures` works with weights", {
  set.seed(258179)
  	weights<-matrix(runif(nrow(seSimCount)*ncol(seSimCount)),nrow=nrow(seSimCount),ncol=ncol(seSimCount))
  ceSimCountW<-ceSimCount
  assay(ceSimCountW,"weights")<-weights
  
	##Test mergeClusters
	expect_silent(ceSimCountW<-makeDendrogram(ceSimCountW))
	expect_message(mergeClusters(ceSimCountW,DEMethod="edgeR"),"Merging will be done on")
	expect_message(mergeClusters(ceSimCountW,DEMethod="edgeR",weights=NULL,forceCalculate=TRUE),"Merging will be done on")
	
	##Test RSEC
	expect_message(RSEC(ceSimCountW,sequential=FALSE,subsample=FALSE,ks=4:15,isCount=TRUE,stopOnError=TRUE),"merging with these parameters did not result in any clusters being merged")
	
	sceCountData<-SingleCellExperiment(simCount)
	reducedDims(sceCountData) <- reducedDims(sceSimDataDimRed)	
	expect_message(RSEC(sceCountData,sequential=TRUE,subsample=FALSE,clusterFunction="pam",k0s=4:15,isCount=TRUE, stopOnError=TRUE),"Merging will be done on")
  
})
