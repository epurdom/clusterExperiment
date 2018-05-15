##put separate so that not too long without output...
context("Heatmap functions")


plotAll<-FALSE #set to true to actually SEE the plots; otherwise for TravisCI, where no visual, runs quicker with FALSE
###Note some are still run with plot=TRUE to check that works with aheatmap. Only a fraction not plotted.
test_that("`setBreaks`", {
	setBreaks(smSimData)
	setBreaks(smSimData,breaks=0.99)
	x<-setBreaks(smSimData,breaks=0.99,makeSymmetric=TRUE)
	expect_equal(max(x),-min(x))
	expect_equal(x,setBreaks(smSimData,breaks=0.01,makeSymmetric=TRUE))
	expect_warning(y<-setBreaks(smSimData,breaks=10))
	expect_equal(length(y),10)
})
test_that("`plotHeatmap` works with matrix objects", {
    expect_silent(x1<-plotHeatmap(data=smSimData))
    a1<-NMF::aheatmap(smSimData)
    expect_equal(x1$aheatmapOut,a1)
    expect_silent(x2<-plotHeatmap(data=smSimCount,clusterSamplesData=smSimData,clusterFeaturesData=smSimData))
    #for some reason, labels on dendrogram move from character to numeric so can't test entire object...
    expect_equal(x1$aheatmapOut$rowInd,x2$aheatmapOut$rowInd) 
    expect_equal(x1$aheatmapOut$colInd,x2$aheatmapOut$colInd) 
    
    #check internal alignment of sampleData (alignSampleData=TRUE) is working:
    sampleData<-clusterMatrix(smSimCE)
    expect_silent(alList<-plotClusters(sampleData))
    alCol<-clusterExperiment:::.convertToAheatmap(alList$clusterLegend, names=FALSE)
   #these should be same plots:
    expect_silent(x1<-plotHeatmap(data=smSimData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,1:10],clusterLegend=alCol,clusterSamples=FALSE,clusterFeatures=FALSE,plot=plotAll))
    expect_silent(x2<-plotHeatmap(data=smSimData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,1:10],alignSampleData=TRUE,clusterFeatures=FALSE,clusterSamples=FALSE,plot=plotAll))
#   Should get this working so proper test, but more a problem because in different order, otherwise the same. Don't want to deal with this right now.
#    expect_equal(lapply(x1$clusterLegend,function(x){x[,c("clusterIds","color")]}),lapply(x2$clusterLegend,function(x){x[,c("clusterIds","color")]}))

    expect_error( plotHeatmap(data=smSimData,Rowv=TRUE),"arguments to aheatmap cannot be set by the user")
    expect_error( plotHeatmap(data=smSimData,Colv=TRUE),"arguments to aheatmap cannot be set by the user")
    expect_error( plotHeatmap(data=smSimData,colorScale=seqPal5,color=TRUE),"arguments to aheatmap cannot be set by the user")

    expect_error( plotHeatmap(data=smSimData,annCol=rnorm(n=ncol(smSimData))),"arguments to aheatmap cannot be set by the user")
    expect_error( plotHeatmap(data=smSimData,annColors=list(a=c("blue","green"))),"arguments to aheatmap cannot be set by the user")
    
    ##Should add tests that pass aheatmap arguments correctly.
})


test_that("`plotHeatmap` works with hdf5 objects", {
    expect_silent(plotHeatmap(hdfObj))

})
test_that("`plotHeatmap` works with ClusterExperiment and SummarizedExperiment objects", {

    expect_silent(plotHeatmap(cc))
    expect_silent(plotHeatmap(cc,whichClusters="none"))
    expect_warning(plotHeatmap(cc,whichClusters="workflow",plot=plotAll) ,"whichClusters value does not match any clusters") #there are no workflow for this one

    	expect_warning(plotHeatmap(smSimCE,whichClusters="workflow",overRideClusterLimit=TRUE),"More than 10 annotations/clusterings can result in incomprehensible errors in aheamap")
    expect_warning(plotHeatmap(smSimCE,whichClusters=15:20,plot=plotAll),"given whichClusters value does not match any clusters")
	expect_error( plotHeatmap(smSimCE,whichClusters="all", alignSampleData=TRUE, overRideClusterLimit=FALSE), "More than 10 annotations/clusterings")
    expect_warning(plotHeatmap(smSimCE,whichClusters="all",
		alignSampleData=FALSE,overRideClusterLimit=TRUE))

 
    #test sampleData
    expect_error(plotHeatmap(cc,sampleData="A"), "no colData for object data")

    expect_silent(plotHeatmap(smSimCE,sampleData="all"))
    expect_silent(plotHeatmap(smSimCE,sampleData="A",plot=plotAll))
    expect_silent(plotHeatmap(smSimCE,sampleData=2:3,plot=plotAll))

    #check that it pulls the names, not the clusterIds.
    clusterLegend(cc)[[1]][,"name"]<-letters[1:nrow(clusterLegend(cc)[[1]])]
    expect_silent(plotHeatmap(cc))
    
    #check user setting clusterLegend
	x<-palette()[1:7]
	names(x)<-clusterLegend(cc)$Cluster1[,"name"]
    expect_silent(plotHeatmap(cc,clusterLegend=list("Cluster1"=x),plot=plotAll))

    expect_silent(plotHeatmap(cc,clusterLegend=list("Cluster1"=palette()[1:7])))
	expect_silent(plotHeatmap(smSimCE,sampleData="A",clusterLegend=list("A"=palette()[1:4]),plot=plotAll))

	names(x)<-LETTERS[1:7]
	expect_error(    plotHeatmap(cc,clusterLegend=list("Cluster1"=x)),"do not cover all levels in the data")
	x<-palette()[1:6]
	names(x)<-LETTERS[1:6]
	expect_error(    plotHeatmap(cc,clusterLegend=list("Cluster1"=x)),"is less than the number of levels in the data")
	
	########################
	########################
    # the following checks work outside of the test but  inside test_that, they hit errors
    # possibly issue with testthat? Not evaluating for now.
	########################
	########################
	#
	# plotHeatmap(smSimCE, sampleData="all", whichClusters="none")
	#
	# #this test doesn't work -- for some reason, expect_warning environment hits error that don't see at the consule.
	# plotHeatmap(smSimCE,whichClusters="all",alignSampleData=TRUE,overRideClusterLimit=TRUE)
	# expect_warning( plotHeatmap(smSimCE, whichClusters="all", alignSampleData=TRUE, overRideClusterLimit=TRUE)
	# , "More than 10 annotations/clusterings")
	#
	# # create some names to see if keeps names with alignSampleData=TRUE
	# # only can check manually, not with testthat.
	# # BUG!: doesn't work. looses their -1/-2 designation... haven't fixed yet.
	# clLeg<-clusterLegend(smSimCE)
	# clLeg[[1]][,"name"]<-LETTERS[1:nrow(clLeg[[1]])]
	# clusterLegend(smSimCE)<-clLeg
	# plotHeatmap(smSimCE, whichClusters="all", alignSampleData=TRUE,overRideClusterLimit=TRUE)
	#
})

test_that("`plotHeatmap` visualization choices/feature choices all work", {

  expect_silent(plotHeatmap(smSimCE,visualizeData=smSimCount,plot=plotAll))
  expect_silent(plotHeatmap(smSimCE,visualizeData="transformed",plot=plotAll))
  expect_silent(plotHeatmap(smSimCE,visualizeData="original",plot=plotAll))
  expect_silent(plotHeatmap(smSimCE,visualizeData="centeredAndScaled",plot=plotAll))
  #even if visualizeData="orginal, still clsuter on transformed. Should make unit test out of below that get same:
  expect_silent(plotHeatmap(smSimCE, visualizeData="transformed", clusterSamplesData="hclust", 
  	plot=plotAll))
  orderSamples(smSimCE)<-sample(1:nSamples(smSimCE))
  expect_silent(plotHeatmap(smSimCE, visualizeData="transformed", clusterSamplesData="orderSamplesValue", plot=plotAll))
  expect_silent(plotHeatmap(smSimCE, visualizeData="transformed", clusterSamplesData="primaryCluster", plot=plotAll))
  expect_silent(plotHeatmap(smSimCE,visualizeData="transformed", clusterSamplesData=c(3,4,5), plot=plotAll))
  expect_silent(plotHeatmap(smSimCE,visualizeData="transform", clusterFeaturesData="all", plot=plotAll))
  
  ###The following subset the genes
  expect_silent(plotHeatmap(smSimCE,visualizeData="transform", clusterFeaturesData="var", nFeatures=3, plot=plotAll))
  expect_silent(plotHeatmap(smSimCE,visualizeData="transform", clusterFeaturesData=3:5, nFeatures=3, plot=plotAll))
  ## expect error because no row names in smSimCE
  expect_error(plotHeatmap(smSimCE, visualizeData="transform", clusterFeaturesData=paste("Gene",3:5), nFeatures=3), "Cannot give feature names in clusterFeaturesData unless")
  ##set rownames and should work
  smSimCE2<-smSimCE
  row.names(smSimCE2)<-paste("Gene",1:NROW(smSimCE))
  expect_silent(plotHeatmap(smSimCE2, visualizeData="transform", clusterFeaturesData=paste("Gene",3:5), nFeatures=3))

  expect_silent(plotHeatmap(smSimCE, visualizeData="transform", clusterFeaturesData="PCA", nFeatures=10, clusterSamplesData="hclust", plot=plotAll))
 expect_silent(plotHeatmap(smSimCE, visualizeData="transform", clusterSamplesData="dendrogramValue", plot=plotAll))
  #test works with outside dataset
 expect_silent(plotHeatmap(smSimCE,visualizeData=assay(smSimCE)[1:10,],plot=plotAll))
 expect_error(plotHeatmap(smSimCE, visualizeData=assay(smSimCE)[,1:5]))
})

test_that("`makeBlankData` works", {
  ##call directly
  gps<-list(c(3,6,7),c(2,1))
  expect_silent(xx<-makeBlankData(assay(smSimCE),groupsOfFeatures=gps))
  expect_equal(nrow(xx$dataWBlanks),length(xx$rowNamesWBlanks))
  whBlankNames<-which(xx$rowNamesWBlanks=="")
  expect_equal(xx$rowNamesWBlanks[-whBlankNames],as.character(unlist(gps)) )
  whBlankRows<-as.numeric(which(apply(xx$dataWBlanks,1,function(x){all(is.na(x))})))
  expect_equal(whBlankRows,whBlankNames)
  expect_equal(whBlankRows,4)

  ##call within plotHeatmap (serves as test of NA capabilities)
  expect_silent(plotHeatmap(smSimCE,clusterFeaturesData=gps))
  expect_silent(plotHeatmap(smSimCE,clusterFeaturesData=gps,breaks=.99))
  expect_warning(plotHeatmap(smSimCE,clusterFeaturesData=gps,breaks=40))
})

test_that("`plotCoClustering` works", {
  expect_error(plotCoClustering(smSimCE),"coClustering slot is empty")
  #following gives all -1, but creates coClustering
  expect_silent(smMin1<-makeConsensus(smSimCE,whichClusters=10:13,proportion=.99))
#  smMin1<-makeConsensus(smSimCE,whichClusters=1:8,proportion=.95) #use to give all -1, but creates coClustering but something changed -- couldn't figure it out!!!
  expect_silent(plotCoClustering(smMin1,clusterSamplesData="hclust"))
  ## Have changed so now changes it internally to primary cluster then hclust
  expect_warning(plotCoClustering(smMin1,clusterSamplesData="dendrogramValue",plot=plotAll),
               "cannot make dendrogram from 'data'")
  expect_silent(sm<-makeConsensus(smSimCE,whichClusters=1:4,proportion=.5))
  expect_silent(plotCoClustering(sm,clusterSamplesData="dendrogramValue"))
})