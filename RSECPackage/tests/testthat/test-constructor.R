context("Constructor")

test_that("`getClusterIndex` works", {
	expect_silent(getClusterIndex(ccSE,whichClusters="workflow"))

	expect_error(getClusterIndex(ccSE,whichClusters="workflow",noMatch="throwError"),"There are no workflow clusterings")
	

	expect_error(getClusterIndex(ccSE,whichClusters="workflow",noMatch="throwError"),"There are no workflow clusterings")
	
	expect_error(getSingleClusterIndex(ccSE,whichCluster="workflow"),"There are no workflow clusterings")
	
	
    expect_message(ccM<-makeConsensus(ceSimCount,proportion=.2,clusterLabel="mySpecialLabel"),"no clusters specified to combine")
    expect_message(ccM<-makeConsensus(ccM,proportion=.2),"no clusters specified to combine")
    expect_silent(indCM<-getClusterIndex(ccM,wh=c("clusterMany")))
    expect_equal(length(indCM),sum(clusterTypes(ccM)=="clusterMany"))
    expect_equal(getClusterIndex(ccM,wh=c("mySpecialLabel","makeConsensus")),c(2,1))
    expect_equal(getClusterIndex(ccM,wh=c("mySpecialLabel","clusterMany")),c(2,indCM))
	expect_equal(getClusterIndex(ccM,wh=c("makeConsensus")),1)
	expect_equal(getClusterIndex(ccM,wh=c("makeConsensus","garbage"),noMatch="silentlyRemove"),1)
})

test_that("whichClusters works with clusterMatrix",{
	 expect_equal(ncol(clusterMatrix(ceSimCount,whichClusters="workflow")),12)
})

test_that("adding clusters work as promised",{

  ## Check coCluster slot with indices
  expect_silent(c5<-addClusterings(ccSE,clusterMatrix(ccSE),clusterTypes="New"))
  expect_silent(coClustering(c5)<-c(4,1))
  expect_silent(c6<-addClusterings(c5,ccSE))
  expect_equal(coClustering(c5),coClustering(c6))
  expect_null(coClustering(addClusterings(ccSE,c5)))
  
})

test_that("removeClusterings work as promised",{

    ## Check coCluster slot with indices
    expect_silent(c4<-addClusterings(ccSE,clusterMatrix(ccSE),clusterTypes="New"))
    expect_error(coClustering(c4)<-c(4,7,1),"CoClustering slot is a vector, but doesn't match indices of clusterMatrix of the object") #can't give more than 4
    expect_silent(coClustering(c4)<-c(4,1))
    ## remove one of the co-clusterings
    expect_warning(c8<-removeClusterings(c4,c(4,3)), "removing clusterings that were used in makeConsensus")
    expect_null(coClustering(c8))
    ## remove separate clustering
    ## remove one of the co-clusterings
    expect_silent(c9<-removeClusterings(c4,c(2)))
    expect_equal(coClustering(c9),c(3,1))
  
})

test_that("assigning unassigned samples works as promised",{
  #also indirectly tests getReducedData!
	expect_silent(assignUnassigned(cc))
	expect_silent(ccUn<-assignUnassigned(cc,whichCluster="Cluster2"))
	expect_true("Cluster2_AllAssigned" %in% clusterLabels(ccUn))
	expect_silent(assignUnassigned(ceSimCount,reduceMethod="b"))
	expect_silent(assignUnassigned(ceSimCount,reduceMethod="mad"))
	expect_silent(ceUn<-assignUnassigned(ceSimCount,reduceMethod="PCA"))
	expect_equal(reducedDimNames(ceUn),"PCA")
	expect_silent(ceUn2<-assignUnassigned(ceSimCount,reduceMethod="PCA",makePrimary=FALSE))
	expect_true(all.equal(primaryCluster(ceSimCount),primaryCluster(ceUn2)))


	#check basic error catching
	cc2<-addClusterings(cc,rep(-1,ncol(cc)),clusterLabel="allUn")
	expect_error(assignUnassigned(cc2,whichCluster="allUn"),"All cells are unassigned, cannot assign them")
	cc2<-addClusterings(cc2,rep(2,ncol(cc2)),clusterLabel="allAss")
	expect_error(assignUnassigned(cc2,whichCluster="allAss"),"No cells are unassigned in the designated cluster")

	#should check whichAssay....
	
})

test_that("workflow functions work",{
  ##########
  #check workflow stuff
  expect_silent(ppC<-addClusterings(cc,cbind(rep(c(-1, 1,2), each=5),rep(c(2, 1,3), each=5)),clusterTypes=c("clusterMany","mergeClusters")))
  expect_equal(dim(workflowClusters(ppC)),c(nSamples(cc),2))
  
  expect_silent(ppC<-addClusterings(cc,cbind(rep(c(-1, 1,2), each=5)),clusterTypes=c("clusterMany")))
  expect_equal(dim(workflowClusters(ppC)),c(nSamples(cc),1))
  
  expect_silent(ppC<-addClusterings(cc,cbind(rep(c(-1, 1,2), each=5),rep(c(2, 3,1), each=5)),clusterTypes=c("clusterMany","mergeClusters.1")))
  expect_equal(dim(workflowClusters(ppC)),c(nSamples(cc),1))
  expect_equal(dim(workflowClusters(ppC,iteration=NA)),c(nSamples(cc),2))
  expect_null(workflowClusters(cc,iteration=NA))
          
  expect_message(ceNew<-makeConsensus(ceSimCount,proportion=0.7))
  expect_message(ceNew<-makeConsensus(ceNew,proportion=0.3,clusterLabel="makeConsensus,v2"))
  expect_equal(clusterLabels(ceNew)[1:2],c("makeConsensus,v2","makeConsensus.1"))
  expect_equal(clusterTypes(ceNew)[1:2],c("makeConsensus","makeConsensus.1"))
  expect_silent(ceNew2<-setToCurrent(ceNew,whichCluster="makeConsensus.1"))
  expect_equal(clusterLabels(ceNew2)[1:2],c("makeConsensus,v2","makeConsensus"))
  expect_equal(clusterTypes(ceNew2)[1:2],c("makeConsensus.2","makeConsensus"))
  expect_silent(ceNew3<-setToCurrent(ceNew2,whichCluster="makeConsensus.2"))
  expect_equal(clusterLabels(ceNew3)[1:2],c("makeConsensus,v2","makeConsensus.3"))
  expect_equal(clusterTypes(ceNew3)[1:2],c("makeConsensus","makeConsensus.3"))

  expect_silent(ceNew4<-setToFinal(ceNew,whichCluster="makeConsensus,v2",clusterLabel="Final Version"))
  expect_equal(primaryClusterIndex(ceNew4),1)
  expect_equal(clusterLabels(ceNew4)[primaryClusterIndex(ceNew4)],"Final Version")
  expect_equal(clusterTypes(ceNew4)[primaryClusterIndex(ceNew4)],"final")
  
  expect_silent(ceNew5<-setToFinal(ceNew,whichCluster="makeConsensus.1",clusterLabel="Final Version"))
  expect_equal(primaryClusterIndex(ceNew5),2)
  expect_equal(clusterLabels(ceNew5)[primaryClusterIndex(ceNew5)],"Final Version")
  expect_equal(clusterTypes(ceNew5)[primaryClusterIndex(ceNew5)],"final")
})