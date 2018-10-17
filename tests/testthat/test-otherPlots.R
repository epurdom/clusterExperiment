context("Non-heatmap related plot functions")



test_that("`plotClusters` works with matrix, ClusterExperiment objects", {

    #test matrix version
    x<-plotClusters(object=clusterMatrix(ceSim))
    expect_equal(dim(clusterMatrix(ceSim)),dim(x$colors))
    expect_equal(dim(clusterMatrix(ceSim)),dim(x$aligned))
    expect_equal(length(x$clusterLegend),ncol(clusterMatrix(ceSim)))

    #test CE version
    x<-plotClusters(ceSim)
    expect_is(x,"ClusterExperiment")
    expect_equal( x,ceSim)

    xx<-plotClusters(ceSim,whichClusters="clusterMany")
    xx2<-plotClusters(object=ceSim,whichClusters="workflow") #only clusterMany values so should be the same
    expect_equal(xx2,xx)

    #check reset -- should add combinations of resetColors and resetNames to make sure works independently.
    par(mfrow=c(1,2)) #so can visually check if desired.
    xx3<-plotClusters(ceSim,resetOrderSamples=TRUE,resetColors=TRUE,resetNames=TRUE)
    plotClusters(xx3,existingColors="all")
    expect_false(isTRUE(all.equal(xx2,xx3))) #not a great test. Doesn't really say whether does it right, just whether it does something!

    nm<-as.numeric(unlist(lapply(clusterLegend(xx3),function(x){x[,"name"]})))
    col<-(unlist(lapply(clusterLegend(xx3),function(x){x[,"color"]})))
    wh<-which(col %in% c("white","grey"))
    expect_equal(match(col[-wh],bigPalette),nm[-wh])
    nmOld<-as.numeric(unlist(lapply(clusterLegend(ceSim),function(x){x[,"name"]})))
    expect_false(isTRUE(all.equal(nm,nmOld)))
    idOld<-as.numeric(unlist(lapply(clusterLegend(ceSim),function(x){x[,"clusterIds"]})))
    idNew<-as.numeric(unlist(lapply(clusterLegend(xx3),function(x){x[,"clusterIds"]})))
    expect_equal(idOld,idNew)

    #check existing colors
    x2<-plotClusters(ceSim,existingColors="all")

    #test -1
    plotClusters(ceSim)

    #CE object with mixture of workflow and other types
    x1<-plotClusters(object=ceSim,whichClusters="workflow",resetColors=TRUE)
    x2<-plotClusters(object=removeClusterings(ceSim,"User"),resetColors=TRUE)
    whP<-getClusterIndex(ceSim,"workflow")
    expect_equal(clusterLegend(x2),clusterLegend(x1)[whP])

    #test specifying indices
    wh<-c(3,4,NCOL(clusterMatrix(ceSim)))
    x3<-plotClusters(ceSim,whichClusters=wh,axisLine=-2,resetColors=TRUE)
    x4<-plotClusters(ceSim,whichClusters=wh[c(3,2,1)],axisLine=-2,resetColors=TRUE)
    expect_false(isTRUE(all.equal(x3,x4)))

    par(mfrow=c(1,1)) #otherwise will affect other tests.
})



test_that("`plotClusters` rerun above tests with colData included", {

  #test matrix version
  x<-plotClusters(object=clusterMatrix(ceSim),colData=as.data.frame(colData(ceSim)))
  expect_equal(ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)),ncol(x$colors))
  expect_equal(ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)),ncol(x$aligned))
  expect_equal(length(x$clusterLegend),ncol(clusterMatrix(ceSim))+ncol(colData(ceSim)))

  #test CE version
  test<- clusterMany(simCount,reduceMethod="PCA",nReducedDims=c(5,10,50), isCount=TRUE,
                     clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE)) #no colData in test
  expect_error(plotClusters(test,colData=as.data.frame(colData(ceSim))),"no colData for object data")
  expect_error(plotClusters(ceSim,colData=as.data.frame(colData(ceSim))),"invalid values for pulling sample data from colData of object")
  plotClusters(ceSim,colData="all")
  par(mfrow=c(1,2))
  x2<-plotClusters(ceSim,colData="all",resetColors=TRUE)
  x1<-plotClusters(ceSim,resetColors=TRUE)


  #check NAs
  naSim<-ceSim
  colData(naSim)[sample(10,1:nrow(naSim)),]<-NA
  plotClusters(naSim,colData=c("A","B"))

  #test the new TRUE option
  plotClusters(naSim,colData=TRUE)

  #this is not working because first one puts -1/-2 last and second puts them first, and so then assigns different colors to the groups
#  expect_equal(x1,x2)
#   par(mfrow=c(1,2))
#   x2<-plotClusters(ceSim,colData="all",resetColors=FALSE)
#   x1<-plotClusters(ceSim,resetColors=FALSE)
  par(mfrow=c(1,1))

})


test_that("plotClustersWorkflow", {
	cc<-clusterMany(mat, ks=c(3,4),nFilterDims=c(10,15),nReducedDims=c(3,4),reduceMethod=c("none","PCA","var"),clusterFunction="pam",
	                       subsample=FALSE, sequential=FALSE,run=TRUE,verbose=FALSE,
	                       isCount=FALSE)
	cc<-makeConsensus(cc,proportion=.7,whichClusters = "clusterMany")
	plotClustersWorkflow(cc)
	plotClustersWorkflow(cc,clusterManyLabels=FALSE)
	plotClustersWorkflow(cc,sortBy="clusterMany")
	plotClustersWorkflow(cc,sortBy="clusterMany",highlightOnTop=FALSE)
	plotClustersWorkflow(cc,highlightOnTop=FALSE)
	plotClustersWorkflow(cc,clusterManyLabels=FALSE,clusterLabels="test")
	expect_error(plotClustersWorkflow(cc,clusterManyLabels=c("1","2"),clusterLabels="test"),"number of cluster labels given in clusterManyLabels")
	expect_error(plotClustersWorkflow(cc,clusterManyLabels=TRUE,clusterLabels=c("A","test")),"number of cluster labels given in clusterLabels")

})


test_that("plotting helpers", {
  convertClusterLegend(smSimCE,output="aheatmap")
  convertClusterLegend(smSimCE,output="plotAndLegend")
  convertClusterLegend(smSimCE,output="matrixNames")
  convertClusterLegend(smSimCE,output="matrixColors")
  convertClusterLegend(smSimCE,output="matrixNames",whichClusters=c("cluster1"))
  convertClusterLegend(smSimCE,output="matrixNames",whichClusters=1:3)
  convertClusterLegend(smSimCE,output="plotAndLegend",whichClusters=c("cluster1"))
  expect_error(convertClusterLegend(smSimCE,output="plotAndLegend",whichClusters=1:3),"given whichClusters indicates more than 1 clustering which is not allowed for option")

  plotClusterLegend(smSimCE)
  plotClusterLegend(smSimCE,whichCluster="cluster1")
    showPalette()
  showPalette(massivePalette)
})



test_that("`plotBarplot` works with matrix, ClusterExperiment objects", {

    #test numeric matrix version
    expect_silent(plotBarplot(object=clusterMatrix(ceSim)[,1:2]))
    #test vector version
    expect_silent(plotBarplot(object=clusterMatrix(ceSim)[,1]))
    #check error
    expect_error(plotBarplot(object=clusterMatrix(ceSim)),"if 'object' a matrix, must contain at most 2 clusters")

    #test CE version with no defaults
    expect_silent(plotBarplot(ceSim))
    #test CE version whichClusters arguments
    expect_silent(plotBarplot(ceSim,whichClusters="workflow"))
    expect_silent(plotBarplot(ceSim,whichClusters="primaryCluster"))
    expect_silent(plotBarplot(ceSim))


    test<-ceSim
    expect_silent(clusterLegend(test)[[1]][,"name"]<-LETTERS[1:nrow(clusterLegend(ceSim)[[1]])])
    #test character matrix version
    expect_silent(plotBarplot(object=convertClusterLegend(test,output="matrixNames")[,1:2]))
    #test character vector version
    expect_silent(plotBarplot(object=convertClusterLegend(test,output="matrixNames")[,1]))
    #test labels argument
    expect_silent(plotBarplot(test,whichClusters=1:2,labels="id"))
    expect_silent(plotBarplot(test,whichClusters=1:2,labels="name"))
    #plotBarplot(ceSim,whichClusters="primaryCluster")

})




test_that("plotReducedDims works",{
	expect_silent(cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reduceMethod="PCA",
	clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
	removeSil=c(TRUE,FALSE)))
	expect_silent(plotReducedDims(cl,legend="bottomright"))
	expect_silent(plotReducedDims(cl,legend=TRUE))
	expect_silent(plotReducedDims(cl,legend=FALSE))
	expect_silent(clusterLegend(cl)[["nReducedDims=10,k=4,findBestK=FALSE,removeSil=TRUE"]][,"name"]<-LETTERS[1:5])
	expect_silent(plotReducedDims(cl,whichCluster="nReducedDims=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE))

	#test on object that doesn't have saved:
	expect_warning(clD<-plotReducedDims(ceSimData,reducedDim="PCA"),"will be run on the FIRST assay")
	expect_equal(NCOL(reducedDim(clD,type="PCA")),2) #default.

	#higher dims.
	expect_silent(plotReducedDims(cl,whichCluster="nReducedDims=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE,whichDims=1:4))
	expect_error(plotReducedDims(cl,whichCluster="nReducedDims=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE,whichDims=158:200),"Invalid value for whichDims: larger than row or column")
	#force it to recalculate:
	expect_warning(plotReducedDims(cl,whichCluster="nReducedDims=10,k=4,findBestK=FALSE,removeSil=TRUE",legend=TRUE,whichDims=51:58),"will be run on the FIRST assay")



})

test_that("plotFeatureBoxplot works",{
	expect_silent(cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reducedDim="PCA",
		clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
		removeSil=c(TRUE,FALSE)))
	expect_silent(clusterLegend(cl)[[1]][,"name"]<-letters[1:nClusters(cl,ignoreUnassigned =FALSE)[1]])
	expect_silent(plotFeatureBoxplot(object=cl,feature=1))
	expect_silent(plotFeatureBoxplot(cc,feature=rownames(cc)[2]))
	expect_silent(plotFeatureBoxplot(cc,plotUnassigned=TRUE,feature=rownames(cc)[2]))
	#check if only 1 non-negative cluster in clustering
	expect_silent(out<-plotFeatureBoxplot(cl[,1:10],whichCluster=2,feature=2,plotUnassigned=FALSE))
	expect_equal(ncol(out$stats),1)
	#check if only 1 cluster in clustering
	expect_silent(out<-plotFeatureBoxplot(cl[,1:20],whichCluster="primary",feature=2))
	expect_equal(ncol(out$stats),1)
})

test_that("plotClustersTable works",{
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
	expect_silent(ceSim<-renameClusters(ceSim,whichCluster=1,val=letters[1:nClusters(ceSim)[1]]))
	expect_silent(ceSim<-subsetByCluster(ceSim,whichCluster=1,c("a","b","d")))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),xlab="Cluster1",margin=2,legend=TRUE))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),xlab=NULL,ylab=NA,margin=0,legend=TRUE))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),xlab=NA,ylab=NULL,margin=1,legend=TRUE))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),margin=NA,ylab="Cluster2",legend=TRUE ))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),margin=NULL,legend=TRUE))
	
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),xlab="Cluster1",margin=2,plotType="bubble"))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),xlab=NULL,ylab=NA,margin=0,plotType="bubble"))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),xlab=NA,ylab=NULL,margin=1,plotType="bubble"))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),margin=NA,ylab="Cluster2",plotType="bubble" ))
	expect_silent(plotClustersTable(ceSim,whichClusters=c(1,2),margin=NULL,plotType="bubble"))
	
	
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

test_that("plotting works with hdf5 assays objects",{
	##plotClusters
    expect_silent(cl1 <- clusterSingle(hdfSCE, reduceMethod="PCA",
            subsample=FALSE, sequential=FALSE,
			mainClusterArgs=list(clusterFunction="pam",clusterArgs=list(k=6)),
			isCount=FALSE))
	expect_silent(plotClusters(cl1))
	
	##plotBarplot
	expect_silent(plotBarplot(cl1))
	
	##plotReducedDims
	expect_silent(plotReducedDims(cl1,legend="bottomright"))

	##plotFeatureBoxplot
	expect_silent(plotFeatureBoxplot(object=cl1,feature=1))

	

})