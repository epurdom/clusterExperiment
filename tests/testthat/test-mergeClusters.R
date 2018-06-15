context("mergeClusters")


test_that("`mergeClusters` works with matrix and ClusterExperiment objects", {
  cl1 <- clusterSingle(smSimData, 
                       subsample=FALSE, sequential=FALSE,
                       mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)), 
											 isCount=FALSE)
  leg<-clusterLegend(cl1)[[primaryClusterIndex(cl1)]]
  leg[,"name"]<-letters[1:6]
  clusterLegend(cl1)[[primaryClusterIndex(cl1)]]<-leg
  clustWithDendro <- makeDendrogram(cl1)
  #matrix version
  mergedList <- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", plotInfo="mergeMethod")
  
	#check plotting types:
  clustMerged <- mergeClusters(clustWithDendro, DEMethod="limma",
	 mergeMethod="none",plotInfo="all")
  clustMerged <- mergeClusters(clustWithDendro, DEMethod="limma", 
		mergeMethod="none", plotInfo="adjP")
  clustMerged <- mergeClusters(clustWithDendro, DEMethod="limma", 
		mergeMethod="none", plotInfo="locfdr")
  expect_warning(clustMerged <- mergeClusters(clustWithDendro, DEMethod="limma",
		 mergeMethod="none", plotInfo="locfdr",showWarnings=TRUE))
  expect_error(clustMerged <- mergeClusters(clustWithDendro, DEMethod="limma",
	 mergeMethod="none", plotInfo="mergeMethod"),
	 "can only plot 'mergeMethod' results if one method is selected")
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", 
		DEMethod="limma", plotInfo="none")

  #check all methods run
  for(method in clusterExperiment:::.availMergeMethods){
	  clustMerged <- mergeClusters(clustWithDendro, DEMethod="limma", 
			mergeMethod=method, plotInfo="mergeMethod")
  }
  
  expect_true("mergeClusters" %in% clusterTypes(clustMerged))
  expect_true("mergeClusters" %in% colnames(clusterMatrix(clustMerged)))

	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", DEMethod="limma",
		 plotInfo="mergeMethod",leafType="samples",plotType="colorblock")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", DEMethod="limma", plotInfo="mergeMethod",leafType="samples",plotType="name")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", DEMethod="limma", plotInfo="mergeMethod",leafType="clusters",plotType="colorblock")
	clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", DEMethod="limma", plotInfo="mergeMethod",leafType="clusters",plotType="name")

  expect_error(mergeClusters(x=transformData(clustWithDendro), DEMethod="limma",
                               cl=primaryCluster(clustWithDendro),plot="none",
                               mergeMethod="adjP",
                               dendro=clustWithDendro@dendro_samples),
                 "Not a valid input dendrogram")

  #test if already exists
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", DEMethod="limma")
  primaryClusterIndex(clustMerged)<-2
  clustMerged<- makeDendrogram(clustMerged)
  clustMerged2<-mergeClusters(clustMerged,mergeMethod="adjP", DEMethod="limma")
  expect_true("mergeClusters.1" %in% clusterTypes(clustMerged2))
  expect_true(!"makeConsensus.1" %in% clusterTypes(clustMerged2))
  expect_true(!"clusterMany.1" %in% clusterTypes(clustMerged2))
  removeClusterings(clustMerged, whichClusters = "mergeClusters")
})


test_that("`mergeClusters` works with HDF5 assay slot",{
    expect_silent(cl1 <- clusterSingle(hdfObj, 
            subsample=FALSE, sequential=FALSE,
			mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
			isCount=FALSE))
    expect_silent(clustWithDendro <- makeDendrogram(cl1))
    expect_message(mergedCE <- mergeClusters(x=clustWithDendro,plot=FALSE,
            mergeMethod="adjP", DEMethod="limma", plotInfo="mergeMethod"),
			"Merging will be done on")
	expect_true(inherits(assay(mergedCE),"DelayedArray"))
								
	
})

test_that("saving merge info works",{
  expect_silent(cl1 <- clusterSingle(smSimData, 
       subsample=FALSE, sequential=FALSE, reduceMethod="none",
	   mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
	   isCount=FALSE))
	   #check cluster is same, otherwise won't get same results...
  expect_equal(primaryCluster(cl1),c(4,1,1,5,1,2,5,1,1,3,3,3,6,6,4,2,2,2))
  #for some reason, clusterSingle not giving reasonable names
  expect_silent(clustWithDendro <- makeDendrogram(cl1,reduceMethod="none"))
  #Dendrogram:
  # --[dendrogram w/ 2 branches and 6 members at h = 1509]
  # |--leaf "3" 
  # `--[dendrogram w/ 2 branches and 5 members at h = 1455]
  # |--[dendrogram w/ 2 branches and 2 members at h = 228]
  # |  |--leaf "1" 
  # |  `--leaf "5" 
  # `--[dendrogram w/ 2 branches and 3 members at h = 337]
  # |--leaf "6" 
  # `--[dendrogram w/ 2 branches and 2 members at h = 259]
  # |--leaf "2" 
  # `--leaf "4"
  #matrix version
  expect_silent(mergedList <- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", plotInfo="none"))
  
  ##check giving nodePropTable
  expect_message(mergedList2<- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1), nodePropTable=mergedList$nodeProp,
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="Storey", plotInfo="none"),"Using existing results of per-node significance")
  expect_equal(mergedList2$nodeProp[,"adjP"],mergedList$nodeProp[,"adjP"])
  expect_silent(mergedListStorey<-mergeClusters(x=transformData(cl1), DEMethod="limma",
                                  cl=primaryCluster(cl1), 
                                  dendro=clustWithDendro@dendro_clusters,
                                  mergeMethod="Storey", plotInfo="mergeMethod"))
  expect_equal(mergedList2$nodeProp[,"Storey"],mergedListStorey$nodeProp[,"Storey"])  
  expect_message(mergedList2Redo<- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1), nodePropTable=mergedList2$nodeProp,
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="none", plotInfo=c("all")),"Using existing results of per-node significance")
  
  #------
  #on ClusterExperiment
  # note need to specify message to see doing it correctly.
  #------
  expect_message(clustMerged <- mergeClusters(clustWithDendro,
	   mergeMethod="adjP",plotInfo="none",plot=FALSE,calculateAll=FALSE),"Merging will be done on")
  expect_equal(clustMerged@merge_dendrocluster_index,clustWithDendro@dendro_index+1)
  expect_equal(clustMerged@merge_dendrocluster_index,clustMerged@dendro_index)
  expect_equal(clustMerged@merge_index,1)
  
  #add to existing with different method 
  expect_message(clustMerged2 <- mergeClusters(clustMerged, mergeMethod="Storey",
  	plotInfo="none",cutoff=0.1,plot=FALSE,calculateAll=FALSE),"Merging will be done on")
  expect_equal(clustMerged2@merge_dendrocluster_index,clustMerged@dendro_index+1)
  expect_equal(clustMerged2@merge_dendrocluster_index,clustMerged2@dendro_index)
  expect_equal(clustMerged2@merge_index,1)
  #rerun previous with different cutoff -- no new calculations
  expect_message(clustMerged3 <- mergeClusters(clustMerged2, mergeMethod="Storey",
  	plotInfo="none",cutoff=0.05,plot=FALSE,calculateAll=FALSE),"Using existing results of per-node significance")
  expect_equal(clustMerged3@merge_dendrocluster_index,clustMerged2@dendro_index+1)
  expect_equal(clustMerged3@merge_dendrocluster_index,clustMerged3@dendro_index)
  expect_equal(clustMerged3@merge_index,1)
  #do it again with higher value that actually merges -- check correct index updates, etc...
  expect_message(clustMerged4 <- mergeClusters(clustMerged3, mergeMethod="Storey",
  	plotInfo="none",cutoff=0.5,plot=FALSE,calculateAll=FALSE),"Using existing results of per-node significance")
  expect_equal(clustMerged4@merge_dendrocluster_index,clustMerged3@dendro_index+1)
  expect_equal(clustMerged4@merge_dendrocluster_index,clustMerged4@dendro_index)
  expect_equal(clustMerged4@merge_index,1)
  expect_equal(clustMerged4@merge_nodeMerge[,"mergeClusterId"],c(NA,NA,1,3,NA))
  expect_equal(clustMerged4@merge_nodeMerge[,"isMerged"],c(FALSE,FALSE,TRUE,TRUE,TRUE))

  #check really gets clusterIds and not names
  expect_silent(clusterLegend(clustMerged3)[["clusterSingle"]][,"name"]<-letters[1:6])
  expect_message(clustMerged5 <- mergeClusters(clustMerged3, mergeMethod="Storey",
  	plotInfo="none",cutoff=0.5,plot=FALSE,calculateAll=FALSE),"Using existing results of per-node significance")
  expect_equal(clustMerged4@merge_nodeMerge,clustMerged5@merge_nodeMerge)
  #helpful for debugging:
  # plotDendrogram(clustMerged,show.node=TRUE,show.tip.label=TRUE)
  # table(clusterMatrix(clustMerged)[,c(1)],clusterMatrix(clustMerged)[,c(2)])
  
  expect_silent(nodeMergeInfo(clustMerged4))
  expect_equal(mergeClusterIndex(clustMerged4),clustMerged4@merge_index)
  expect_equal(mergeCutoff(clustMerged4),0.5)
  expect_equal(mergeMethod(clustMerged4),"Storey")
  
  #check if can calculate all, but do nothing else
  expect_message(clustMergedAll<-mergeClusters(clustWithDendro, mergeMethod="none",plotInfo="none",cutoff=0.5,plot=FALSE,calculateAll=TRUE),"Merging will be done on")
  expect_false(is.na(clustMergedAll@merge_dendrocluster_index))
  expect_true(is.na(clustMergedAll@merge_index))
  expect_false(is.null(nodeMergeInfo(clustMergedAll)))
  
  #should erase merge info if call dendrogram
  expect_silent(clustMergedErase<-makeDendrogram(clustMerged5))
  expect_true(is.na(clustMergedErase@merge_index))
  expect_true(is.na(clustMergedErase@merge_dendrocluster_index))
  
  #should erase merge info if call dendrogram
  expect_silent(clustMergedErase2<-makeDendrogram(clustMergedAll))
  expect_true(is.na(clustMergedErase2@merge_index))
  expect_true(is.na(clustMergedErase2@merge_dendrocluster_index))

  #test getMergeCorresp
  #node clustMerge
  mgCl<-clusterMatrix(clustMerged)[,clustMerged@merge_index]
  ogCl<-clusterMatrix(clustMerged)[,clustMerged@merge_dendrocluster_index]
  expect_silent(mc<-getMergeCorrespond(clustMerged))
  expect_equal(length(mc),length(unique(mgCl[mgCl>0])))
  expect_silent(mc<-getMergeCorrespond(clustMerged,by="original"))
  expect_equal(length(mc),length(unique(ogCl[ogCl>0])))
  
  expect_error(getMergeCorrespond(clustMergedAll),"there is no merge clustering in this object")
  #check saves even if just plotInfo
  clustMerged <- mergeClusters(clustWithDendro, mergeMethod="none",plotInfo="all")
  expect_false(is.null(clustMerged@merge_nodeProp))
})


test_that("logFC works",{
  expect_silent(cl1 <- clusterSingle(smSimData, 
       subsample=FALSE, sequential=FALSE, reduceMethod="none",
	   mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
	   isCount=FALSE))
	   #check cluster is same, otherwise won't get same results...
  expect_silent(clustWithDendro <- makeDendrogram(cl1,reduceMethod="none"))
  #--------
  #matrix version 
  #--------
  #only FC
  expect_silent(mergedList <- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", cutoff=0.1, logFCcutoff=1,calculateAll=FALSE,plotInfo="none",plot=FALSE))
  expect_true(all(c(clusterExperiment:::.availMergeMethods,"adjP_1.0") %in% colnames(mergedList$nodeProp)))

  #calculate All 
  expect_silent(mergedListAll <- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", cutoff=0.1, logFCcutoff=1,calculateAll=TRUE,plotInfo="none",plot=FALSE))
  expect_true(all(c(clusterExperiment:::.availMergeMethods,"adjP_1.0") %in% colnames(mergedListAll$nodeProp)))
  
  #check calculate all does FC it even if not given as the mergeMethod
  # (try both "none" and "Storey")
  expect_silent(mergedListAll2 <- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="none",
							   logFCcutoff=1,plotInfo="none",plot=FALSE))
  expect_true(all(c(clusterExperiment:::.availMergeMethods,"adjP_1.0") %in% colnames(mergedListAll2$nodeProp)))
  expect_silent(mergedListAll2 <- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1),
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="Storey",
							   logFCcutoff=1,plotInfo="none",plot=FALSE))
  expect_true(all(c(clusterExperiment:::.availMergeMethods,"adjP_1.0") %in% colnames(mergedListAll2$nodeProp)))
  
  ##check giving nodePropTable but need new calculation ("Storey")
  expect_silent(mergedList2<- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1), nodePropTable=mergedList$nodeProp,
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="Storey", plotInfo="none",plot=FALSE))
  expect_equal(mergedList2$nodeProp[,"adjP_1.0"],mergedList$nodeProp[,"adjP_1.0"])

  ##check giving nodePropTable but need new calculation (newFC)
  ##Also check that actually gives different merging (with logFC10 should)
  expect_silent(mergedList2<- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1), nodePropTable=mergedList$nodeProp,
                              dendro=clustWithDendro@dendro_clusters, cutoff=0.1,
                              mergeMethod="adj", logFCcutoff=10,
							  plotInfo="none",plot=FALSE))
  expect_equal(mergedList2$nodeProp[,"adjP_1.0"],mergedList$nodeProp[,"adjP_1.0"])
  expect_true("adjP_10.0" %in% names(mergedList2$nodeProp))
  expect_false(all(mergedList2$nodeProp[,"adjP_10.0"] == mergedList$nodeProp[,"adjP"]))
  expect_false(all(mergedList2$nodeMerge[,"isMerged"] == mergedList$nodeMerge[,"isMerged"]))
  expect_equal(length(unique(mergedList2$clustering)),3)
  expect_equal(length(unique(mergedList$clustering)),5)

  ##check giving nodePropTable but DON'T need new calculation ("Storey")
  expect_message(mergedList3<- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1), nodePropTable=mergedListAll$nodeProp,
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="Storey", plotInfo="none",plot=FALSE),"Using existing results of per-node significance")
  expect_equal(mergedList2$nodeProp[,"adjP_1.0"],mergedList$nodeProp[,"adjP_1.0"])

  ##check giving nodePropTable but DON'T need new calculation (logFC)
  expect_message(mergedList2<- mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1), nodePropTable=mergedListAll$nodeProp,
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", logFCcutoff=1, plotInfo="none",plot=FALSE),"Using existing results of per-node significance")
  expect_equal(mergedList2$nodeProp[,"adjP_1.0"],mergedList$nodeProp[,"adjP_1.0"])
  
  ###Test plotting on matrix version:
  expect_silent(mergeClusters(x=transformData(cl1), DEMethod="limma",
                              cl=primaryCluster(cl1), 
                              dendro=clustWithDendro@dendro_clusters,
                              mergeMethod="adjP", logFCcutoff=10, plotInfo="adjP_10.0"))
  
  
  #----------
  #on ClusterExperiment
  #----------
  expect_message(clustMerged <- mergeClusters(clustWithDendro,
	   mergeMethod="adjP", plotInfo="none", plot=FALSE, 
	   calculateAll=FALSE,cutoff=0.1,logFCcutoff=1),"Merging will be done on")
  expect_true(all(c("adjP_1.0") %in% colnames(clustMerged@merge_nodeProp)))
  
  #add to existing with different method
  expect_message(clustMerged2 <- mergeClusters(clustMerged, mergeMethod="Storey",
  	plotInfo="none",cutoff=0.1,plot=FALSE,calculateAll=FALSE),"Merging will be done on")
  expect_true(all(c("adjP_1.0") %in% colnames(clustMerged2@merge_nodeProp)))
  
  #rerun previous with different cutoff but same logFC cutoff -- no new calculations needed
  expect_message(clustMerged3 <- mergeClusters(clustMerged2, mergeMethod="adjP", logFCcutoff=1,plotInfo="none",cutoff=0.05,plot=FALSE,calculateAll=FALSE),"Using existing results of per-node significance")
  expect_true(all(c("adjP_1.0") %in% colnames(clustMerged3@merge_nodeProp)))

  #redo with different logFC --- pick logFC large enought that changes the proportions to make sure actually merge on right value
  expect_message(clustMerged4 <- mergeClusters(clustMerged, mergeMethod="adjP", 
                                               logFCcutoff=10,plotInfo="none",cutoff=0.1,plot=FALSE,
                                               calculateAll=FALSE,clusterLabel="merge, logFC"),
                 "Merging will be done on")
  expect_false(all(clustMerged4@merge_nodeProp[,"adjP_10.0"] == clustMerged4@merge_nodeProp[,"adjP"]))
  expect_false(all(clustMerged4@merge_nodeMerge[,"isMerged"] == clustMerged@merge_nodeMerge[,"isMerged"]))
  expect_equal(length(unique(primaryCluster(clustMerged4))),3)
  expect_equal(length(unique(primaryCluster(clustMerged))),5)

  ##Check labels not getting corrupted when redo
  expect_equal(primaryClusterLabel(clustMerged4),"merge, logFC")
  expect_message(clustMerged5 <- mergeClusters(clustMerged4, mergeMethod="adjP", logFCcutoff=10,plotInfo="none",
                                               cutoff=0.1,plot=FALSE,calculateAll=FALSE,clusterLabel="merge, logFC"),"Merging will be done on")
  expect_equal(primaryClusterLabel(clustMerged5),"merge, logFC")
  
  
  ###Test plotting on CE version:
  expect_message(clustMerged10 <- mergeClusters(clustWithDendro,
	   mergeMethod="adjP", plot=TRUE,plotInfo="adjP_10.0", 
	   calculateAll=FALSE,cutoff=0.1,logFCcutoff=10),"Merging will be done on")

  expect_silent(plotDendrogram(clustMerged10,mergeInfo="adjP_10.0"))
})


test_that("`mergeClusters` preserves the colData and rowData of SE", {

  expect_silent(cl <- clusterSingle(smSimSE, 
       subsample=FALSE, sequential=FALSE,                      
	   mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
	   isCount=FALSE))
  expect_silent(clD <- makeDendrogram(cl))
  expect_message(cl <- mergeClusters(clD, mergeMethod = "adjP"))
  expect_equal(colData(cl),colData(smSimSE))
  expect_equal(rownames(cl),rownames(smSimSE))
  expect_equal(colnames(cl),colnames(smSimSE))
  expect_equal(metadata(cl),metadata(smSimSE))
  expect_equal(rowData(cl),rowData(clD))

})


test_that("`mergeClusters` works with unassignedSamples", {

  expect_silent(clustWithDendro <- makeDendrogram(ceSim,unassignedSamples = c("outgroup")))

	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP",
		 plotInfo="mergeMethod",leafType="samples",plotType="colorblock"))
	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP",
		 plotInfo="mergeMethod",leafType="samples",plotType="name"))
	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP",
		 plotInfo="mergeMethod",leafType="clusters",plotType="colorblock"))
	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",plotType="name"))

  expect_silent(clustWithDendro <- makeDendrogram(ceSim,reduceMethod="mad",unassignedSamples = c("cluster")))

	expect_warning(mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",plotType="colorblock"),"You cannot set 'leafType' to 'samples' in plotting mergeClusters unless the dendrogram was made with unassigned/missing")
	expect_warning(mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="samples",plotType="name"),"You cannot set 'leafType' to 'samples' in plotting mergeClusters unless the dendrogram was made with unassigned/missing")
	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",plotType="colorblock","You cannot set 'leafType' to 'samples' in plotting mergeClusters unless the dendrogram was made with unassigned/missing"))
	expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP", plotInfo="mergeMethod",leafType="clusters",plotType="name","You cannot set 'leafType' to 'samples' in plotting mergeClusters unless the dendrogram was made with unassigned/missing"))


})

test_that("cluster labels not being internally changed from user input",{
  expect_silent(clustWithDendro <- makeDendrogram(ceSim,unassignedSamples = c("outgroup")))
  
  expect_message(clustMerged <- mergeClusters(clustWithDendro, mergeMethod="adjP",
                                              plotInfo="mergeMethod",leafType="samples",plotType="colorblock"))
  
  
  })