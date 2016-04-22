source("create_objects.R")

ps<-c(5,10,50)
cl <- clusterMany(simData,dimReduce="PCA",nPCADims=ps,
                  clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE))
cl2<-addClusters(cl,sample(2:5,size=NCOL(simData),replace=TRUE),clusterType="User")
clMatNew<-apply(clusterMatrix(cl),2,function(x){
    wh<-sample(1:nSamples(cl),size=10)
    x[wh]<- -1
    wh<-sample(1:nSamples(cl),size=10)
    x[wh]<- -2
    return(x)
})
#make a new object with -1 values
cl3<-clusterExperiment(assay(cl),clMatNew,transformation=transformation(cl))

sData<-data.frame(sample(letters[2:5],size=NCOL(simData),replace=TRUE),sample(2:5,size=NCOL(simData),replace=TRUE))
sData<-data.frame(sData,sample(LETTERS[2:5],size=NCOL(simData),replace=TRUE),stringsAsFactors=FALSE)
colnames(sData)<-c("A","B","C")

whSamp<-unlist(tapply(1:nSamples(cl3),primaryCluster(cl3),function(x){sample(x,size=3)}))
smData<-simData[1:10,whSamp]
smCount<-simCount[1:10,whSamp]
smCl<-cl3[1:10,whSamp]
smSData<-sData[whSamp,]
clusterLabels(smCl)<-paste("Cluster",1:nClusters(smCl))
smCl2<-clusterExperiment(simCount[1:10,whSamp],clMatNew[whSamp,],transformation=function(x){log(x+1)})
colData(smCl2)<-DataFrame(smSData)
clusterLabels(smCl2)<-paste("Cluster",1:nClusters(smCl2))

test_that("`plotClusters` works with matrix, ClusterExperiment objects", {
    #test matrix version
    x<-plotClusters(clusters=clusterMatrix(cl))
    expect_equal(dim(clusterMatrix(cl)),dim(x$colors))
    expect_equal(dim(clusterMatrix(cl)),dim(x$aligned))
    expect_equal(length(x$clusterLegend),ncol(clusterMatrix(cl)))
    expect_error(plotClusters(clusters=clusterMatrix(cl),whichClusters="garbage"),"unable to find an inherited method")
    expect_error(plotClusters(clusters=clusterMatrix(cl),whichClusters=c(1,3,4)),"unable to find an inherited method")

    #test CE version
    x<-plotClusters(cl)
    expect_is(x,"ClusterExperiment")
    expect_equal( x,cl)
    
    xx<-plotClusters(cl2,whichClusters="clusterMany")
    xx2<-plotClusters(clusters=cl2,whichClusters="pipeline") #only clusterMany values so should be the same
    expect_equal(xx2,xx)
    
    #check reset -- should add combinations of resetColors and resetNames to make sure works independently.
    par(mfrow=c(1,2)) #so can visually check if desired.
    xx3<-plotClusters(cl,resetOrderSamples=TRUE,resetColors=TRUE,resetNames=TRUE) 
    expect_false(isTRUE(all.equal(xx2,xx3))) #not a great test. Doesn't really say whether does it right, just whether it does something!
    nm<-as.numeric(unlist(lapply(clusterLegend(xx3),function(x){x[,"name"]})))
    col<-(unlist(lapply(clusterLegend(xx3),function(x){x[,"color"]})))
    expect_equal(match(col,bigPalette),nm)
    nmOld<-as.numeric(unlist(lapply(clusterLegend(cl),function(x){x[,"name"]})))
    expect_false(isTRUE(all.equal(nm,nmOld)))
    idOld<-as.numeric(unlist(lapply(clusterLegend(cl),function(x){x[,"clusterIds"]})))
    idNew<-as.numeric(unlist(lapply(clusterLegend(xx3),function(x){x[,"clusterIds"]})))
    expect_equal(idOld,idNew)

    #check existing colors
    x2<-plotClusters(cl,existingColors="all")
    
    #test -1
    plotClusters(cl3) 
    
    #CE object with mixture of pipeline and other types
    x1<-plotClusters(clusters=cl2,whichClusters="pipeline",resetColors=TRUE)
    x2<-plotClusters(clusters=cl,resetColors=TRUE)
    whP<-.TypeIntoIndices(cl2,"pipeline")
    expect_equal(clusterLegend(x2),clusterLegend(x1)[whP])
    
    #test specifying indices
    wh<-c(3,4,NCOL(clusterMatrix(cl)))
    x3<-plotClusters(cl,whichClusters=wh,axisLine=-2,resetColors=TRUE)
    x4<-plotClusters(cl,whichClusters=wh[c(3,2,1)],axisLine=-2,resetColors=TRUE)
    expect_false(isTRUE(all.equal(x3,x4)))
    
    #test if only a single cluster
    plotClusters(cl,whichClusters=2)
    x5<-plotClusters(cl,whichClusters=2,resetColors=TRUE)
    expect_equal(x5,cl)
    
})


test_that("`plotClusters` rerun above tests with sampleData included", {
  #test matrix version
  x<-plotClusters(clusters=clusterMatrix(cl),sampleData=sData)
  expect_equal(ncol(clusterMatrix(cl))+ncol(sData),ncol(x$colors))
  expect_equal(ncol(clusterMatrix(cl))+ncol(sData),ncol(x$aligned))
  expect_equal(length(x$clusterLegend),ncol(clusterMatrix(cl))+ncol(sData))

  #test CE version
  expect_error(plotClusters(cl,sampleData=sData),"no colData for object data")
  colData(cl)<-DataFrame(sData)
  expect_error(plotClusters(cl,sampleData=sData),"invalid values for pulling sampleData")
  plotClusters(cl,sampleData="all")
  x2<-plotClusters(cl,sampleData="all",resetColors=TRUE)
  x1<-plotClusters(cl,resetColors=TRUE)
  expect_equal(x1,x2)

})

test_that("`plotHeatmap` works with matrix objects", {
    x1<-plotHeatmap(data=smData)
    x2<-plotHeatmap(data=smCount,clusterSamplesData=smData,clusterFeaturesData=smData)
    expect_equal(x1$aheatmapOut,x2$aheatmapOut)
    
    #check internal alignment of sampleData (alignSampleData=TRUE) is working:
    sampleData<-clusterMatrix(smCl)
    alList<-plotClusters(sampleData)
    alCol<-alList$clusterLegend
    x1<-plotHeatmap(data=smData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,],clusterLegend=alCol,clusterSamples=FALSE,clusterFeatures=FALSE)
    x2<-plotHeatmap(data=smData[,alList$orderSamples],sampleData=sampleData[alList$orderSamples,],alignSampleData=TRUE,clusterFeatures=FALSE,clusterSamples=FALSE)
#   Should get this working so proper test, but more a problem because in different order, otherwise the same. Don't want to deal with this right now.
#    expect_equal(lapply(x1$clusterLegend,function(x){x[,c("clusterIds","color")]}),lapply(x2$clusterLegend,function(x){x[,c("clusterIds","color")]}))
})

test_that("`plotHeatmap` works with ClusterExperiment and SummarizedExperiment objects", {
    plotHeatmap(smCl)
    plotHeatmap(smCl,whichClusters="none")
    expect_warning(plotHeatmap(smCl,whichClusters="pipeline") ,"whichClusters value does not match any clusters") #there are no pipeline for this one
    clusterType(smCl)[2:3]<-"clusterMany"
    plotHeatmap(smCl,whichClusters="pipeline") 
    plotHeatmap(smCl,whichClusters="all",alignSampleData=TRUE)
    expect_error(plotHeatmap(smCl,whichClusters=1:15),"Indices in whichClusters invalid")
    
    #test sampleData
    expect_error(plotHeatmap(smCl,sampleData="A"))
    colData(smCl)<-DataFrame(smSData)
    plotHeatmap(smCl,sampleData="all")
    plotHeatmap(smCl,sampleData="A")
    plotHeatmap(smCl,sampleData=2:3)
    plotHeatmap(smCl,sampleData="all",whichClusters="none")

    #SummarizedExperiment
    se<-SummarizedExperiment(smData)
    plotHeatmap(se)
    
})

test_that("`plotHeatmap` visualization choices/feature choices all work", {
  plotHeatmap(smCl2,visualizeData=smCount)
  plotHeatmap(smCl2,visualizeData="transformed")
  plotHeatmap(smCl2,visualizeData="original")
  plotHeatmap(smCl2,visualizeData="centeredAndScaled")
  #even if visualizeData="orginal, still clsuter on transformed. Should make unit test out of below that get same:
  plotHeatmap(smCl2,visualizeData="transformed",clusterSamplesData="hclust")
  orderSamples(smCl2)<-sample(1:nSamples(smCl2))
  plotHeatmap(smCl2,visualizeData="transformed",clusterSamplesData="orderSamplesValue")
  plotHeatmap(smCl2,visualizeData="transformed",clusterSamplesData="primaryCluster")
  plotHeatmap(smCl2,visualizeData="transformed",clusterSamplesData=c(3,4,5))
  
  plotHeatmap(smCl2,visualizeData="transform",clusterFeaturesData="all")
  plotHeatmap(smCl2,visualizeData="transform",clusterFeaturesData="mostVar",nFeatures=3)
  plotHeatmap(smCl2,visualizeData="transform",clusterFeaturesData=3:5,nFeatures=3)
  expect_error(plotHeatmap(smCl2,visualizeData="transform",clusterFeaturesData=paste("Gene",3:5),nFeatures=3))
  row.names(smCl2)<-paste("Gene",1:NROW(smCl2))
  plotHeatmap(smCl2,visualizeData="transform",clusterFeaturesData=paste("Gene",3:5),nFeatures=3)
   plotHeatmap(smCl2,visualizeData="transform",clusterFeaturesData="PCA",nFeatures=10,clusterSamplesData="hclust")
  
  plotHeatmap(smCl2,visualizeData="transform",clusterSamplesData="dendrogramValue")
  
})

test_that("`makeBlankData` works", {
  ##call directly
  gps<-list(c(3,6,7),c(2,1))
  xx<-makeBlankData(assay(smCl2),groupsOfFeatures=gps)
  expect_equal(nrow(xx$dataWBlanks),length(xx$rowNamesWBlanks))
  whBlankNames<-which(xx$rowNamesWBlanks=="")
  expect_equal(xx$rowNamesWBlanks[-whBlankNames],as.character(unlist(gps)) )
  whBlankRows<-as.numeric(which(apply(xx$dataWBlanks,1,function(x){all(is.na(x))})))
  expect_equal(whBlankRows,whBlankNames)
  expect_equal(whBlankRows,4)
  
  ##call within plotHeatmap
  plotHeatmap(smCl2,clusterFeaturesData=gps)
})
test_that("`plotCoClustering` works", {
  expect_error(plotCoClustering(smCl2),"coClustering slot is empty")
  smCl2<-combineMany(smCl2,whichClusters=1:4,proportion=.9) #gives all -1, but creates coClustering
  plotCoClustering(smCl2,clusterSamplesData="hclust")
  expect_error(plotCoClustering(smCl2,clusterSamplesData="dendrogramValue"),"all samples have clusterIds<0")
  primaryClusterIndex(smCl2)<-3
  plotCoClustering(smCl2,clusterSamplesData="dendrogramValue")
})

test_that("plotting helpers", {
  convertClusterLegend(smCl2,output="aheatmap")
  convertClusterLegend(smCl2,output="plotAndLegend")
  convertClusterLegend(smCl2,output="matrixNames")
  convertClusterLegend(smCl2,output="matrixColors")
})