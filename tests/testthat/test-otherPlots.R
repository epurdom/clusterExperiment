context("Plotting Functions")


test_that("plotting helpers", {
  expect_silent(convertClusterLegend(smSimCE,output="aheatmap"))
  expect_silent(convertClusterLegend(smSimCE,output="plotAndLegend"))
  expect_silent(convertClusterLegend(smSimCE,output="matrixNames"))
  expect_silent(convertClusterLegend(smSimCE,output="matrixColors"))
  expect_silent(convertClusterLegend(smSimCE,output="matrixNames",whichClusters=c("cluster1")))
  expect_silent(convertClusterLegend(smSimCE,output="matrixNames",whichClusters=1:3))
  expect_silent(convertClusterLegend(smSimCE,output="plotAndLegend",whichClusters=c("cluster1")))
  expect_error(convertClusterLegend(smSimCE,output="plotAndLegend",whichClusters=1:3),"given whichClusters indicates more than 1 clustering which is not allowed for option")

  expect_silent(plotClusterLegend(smSimCE))
  expect_silent(plotClusterLegend(smSimCE,whichCluster="cluster1"))
  expect_silent(showPalette())
  expect_silent(showPalette(massivePalette))
})

test_that("`plotBarplot` works with matrix, ClusterExperiment objects", {

    #test numeric matrix version
    expect_silent(plotBarplot(object=clusterMatrix(ceSimCount)[,1:2]))
    #test vector version
    expect_silent(plotBarplot(object=clusterMatrix(ceSimCount)[,1]))
    #check error
    expect_error(plotBarplot(object=clusterMatrix(ceSimCount)),"if 'object' a matrix, must contain at most 2 clusters")

    #test CE version with no defaults
    expect_silent(plotBarplot(ceSimCount))
    #test CE version whichClusters arguments
    expect_silent(plotBarplot(ceSimCount,whichClusters="workflow"))
    expect_silent(plotBarplot(ceSimCount,whichClusters="primaryCluster"))
    expect_silent(plotBarplot(ceSimCount))


    test<-ceSimCount
    expect_silent(clusterLegend(test)[[1]][,"name"]<-LETTERS[1:nrow(clusterLegend(ceSimCount)[[1]])])
    #test character matrix version
    expect_silent(plotBarplot(object=convertClusterLegend(test,output="matrixNames")[,1:2]))
    #test character vector version
    expect_silent(plotBarplot(object=convertClusterLegend(test,output="matrixNames")[,1]))
    #test labels argument
    expect_silent(plotBarplot(test,whichClusters=1:2,labels="id"))
    expect_silent(plotBarplot(test,whichClusters=1:2,labels="name"))
    #plotBarplot(ceSimCount,whichClusters="primaryCluster")

})

test_that("plotReducedDims works",{
	expect_silent(cl <- clusterMany(simData, 
        nReducedDims=c(5, 10, 50), reduceMethod="PCA",
	    clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
	    removeSil=c(TRUE,FALSE), 
        verbose=FALSE,makeMissingDiss=TRUE))
	expect_silent(plotReducedDims(cl,legend="bottomright"))
	expect_silent(plotReducedDims(cl,legend=TRUE))
	expect_silent(plotReducedDims(cl,legend=FALSE))
    whCl<-grep("k=4,nReducedDims=10,findBestK=FALSE,removeSil=TRUE",
        clusterLabels(cl))
    expect_true(length(whCl)==1)
    whCl<-clusterLabels(cl)[whCl] #so is character value, to test.
	expect_silent(clusterLegend(cl)[[whCl]][,"name"]<-LETTERS[1:5])
	expect_silent(plotReducedDims(cl, whichCluster = whCl, legend=TRUE))

	#test on object that doesn't have saved:
	expect_warning(clD<-plotReducedDims(ceSimData,reducedDim="PCA"),
        "will be run on the FIRST assay")
	expect_equal(NCOL(reducedDim(clD,type="PCA")),2) #default.

	#higher dims.
	expect_silent(plotReducedDims(cl,whichCluster=whCl,
        legend=TRUE,whichDims=1:4))
	expect_error(plotReducedDims(cl,whichCluster=whCl,
        legend=TRUE,whichDims=158:200),
        "Invalid value for whichDims: larger than row or column")
	#force it to recalculate:
	expect_warning(plotReducedDims(cl,whichCluster=whCl,legend=TRUE,
        whichDims=51:58),"will be run on the FIRST assay")



})

test_that("plotFeatureBoxplot works",{
	expect_silent(cl <- clusterMany(simData, 
        nReducedDims=c(5, 10, 50), reducedDim="PCA",
		clusterFunction="pam", ks=2:4, 
        findBestK=c(TRUE,FALSE), verbose=FALSE,
		removeSil=c(TRUE,FALSE), makeMissingDiss=TRUE))
	expect_silent(clusterLegend(cl)[[1]][,"name"] <- 
        letters[1:nClusters(cl,ignoreUnassigned =FALSE)[1]])
	expect_silent(plotFeatureBoxplot(object=cl,feature=1))
	expect_silent(plotFeatureBoxplot(cc,feature=rownames(cc)[2]))
	expect_silent(plotFeatureBoxplot(cc,
        plotUnassigned=TRUE,feature=rownames(cc)[2]))
	#check if only 1 non-negative cluster in clustering
	expect_silent(out<-plotFeatureBoxplot(cl[,1:10],
        whichCluster=2,feature=2,plotUnassigned=FALSE))
	expect_equal(ncol(out$stats),1)
	#check if only 1 cluster in clustering
	expect_silent(out<-plotFeatureBoxplot(cl[,1:20],
        whichCluster="primary",feature=2))
	expect_equal(ncol(out$stats),1)
})

test_that("plotClustersTable works",{
    skip_on_os("windows")
	#test where should be diagonal
	expect_silent(plotClustersTable(cc,whichClusters=c(1,2)))
	expect_silent(plotClustersTable(cc,whichClusters=c(1,2),ignoreUnassigned=TRUE,margin=2))
	expect_silent(plotClustersTable(cc,whichClusters=c(1,2),ignoreUnassigned=TRUE,margin=0))
	expect_silent(plotClustersTable(cc,whichClusters=c(1,2),ignoreUnassigned=TRUE,margin=NA))
	expect_silent(plotClustersTable(tableClusters(cc,whichClusters=c(1,2))))

	
	#the following gives a wicked output which ignoreUnassigned=TRUE: zero overlap because some in one cluster are all -1 in makeConsensus so NaN value in proportion and only single cluster in makeConsensus
	expect_silent(cc<-clusterMany(mat, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),reduceMethod=c("none","PCA","var"),clusterFunction="pam",
	                       subsample=FALSE, sequential=FALSE,run=TRUE,verbose=FALSE,
	                       isCount=FALSE))
	expect_message(cc<-makeConsensus(cc,proportion=0.7),"no clusters specified to combine")
	expect_error(plotClustersTable(cc,whichClusters=c(1,2),ignoreUnassigned=TRUE) ,"Cannot create heatmap when there is only 1 column or row in the table")
	expect_silent(plotClustersTable(cc,whichClusters=c(1,2),ignoreUnassigned=TRUE,plotType="bubble")) #gives a 
	
	
	#test more complicated
	#so different numbers of clusters in the two clusters
	expect_silent(ceSimCount<-renameClusters(ceSimCount,whichCluster=1,val=letters[1:nClusters(ceSimCount)[1]]))
	expect_silent(ceSimCount<-subsetByCluster(ceSimCount,whichCluster=1,c("a","b","d")))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),xlab="Cluster1",margin=2,legend=TRUE))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),xlab=NULL,ylab=NA,margin=0,legend=TRUE))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),xlab=NA,ylab=NULL,margin=1,legend=TRUE))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),margin=NA,ylab="Cluster2",legend=TRUE ))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),margin=NULL,legend=TRUE))
	
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),xlab="Cluster1",margin=2,plotType="bubble"))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),xlab=NULL,ylab=NA,margin=0,plotType="bubble"))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),xlab=NA,ylab=NULL,margin=1,plotType="bubble"))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),margin=NA,ylab="Cluster2",plotType="bubble" ))
	expect_silent(plotClustersTable(ceSimCount,whichClusters=c(1,2),margin=NULL,plotType="bubble"))
	
	
})

test_that("plotFeatureScatter works",{
	expect_silent(plotFeatureScatter(object=cc,features=c(1,2),whichCluster=1,pch=19))
	expect_silent(plotFeatureScatter(object=cc,features=c(1,2,3),whichCluster=1,pch=19))
	expect_error(plotFeatureScatter(object=cc,features=c("Gene1","Gene4"),whichCluster=1),"not all of features match one")

	expect_silent(plotFeatureScatter(object=cc,features=c("Gene 1","Gene 4"),whichCluster=1,pch=19))
	expect_silent(plotFeatureScatter(object=cc,features=c("Gene 1","Gene 4","Gene 10"),whichCluster=1,pch=19))

	expect_silent(plotFeatureScatter(object=cc,features=c("Gene 1","Gene 4","Gene 10"),whichCluster=1,pch=19,plotUnassigned=FALSE))
	expect_silent(plotFeatureScatter(object=cc,features=c("Gene 1","Gene 4"),unassignedColor="black",whichCluster=1,pch=19,legend="topright"))

	cc2<-cc
	rownames(cc2)<-NULL
	plotFeatureScatter(object=cc2,features=c(1,2),whichCluster=1,pch=19)
	plotFeatureScatter(object=cc2,features=c(1,2,3),whichCluster=1,pch=19)

})

