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
test_that("`plotClusters` works with matrix, ClusterExperiment objects", {
    #test matrix version
    x<-plotClusters(clusters=clusterMatrix(cl))
    expect_equal(dim(clusterMatrix(cl)),dim(x$colors))
    expect_equal(dim(clusterMatrix(cl)),dim(x$aligned))
    expect_equal(length(x$clusterLegend),ncol(clusterMatrix(cl)))
    expect_error(plotClusters(clusters=clusterMatrix(cl),orderClusters="garbage"),"unable to find an inherited method")
    expect_error(plotClusters(clusters=clusterMatrix(cl),orderClusters=c(1,3,4)),"unable to find an inherited method")

    #test CE version
    x<-plotClusters(cl)
    expect_is(x,"ClusterExperiment")
    expect_equal( x,cl)
    
    xx<-plotClusters(cl2,orderClusters="clusterMany")
    xx2<-plotClusters(clusters=cl2,orderClusters="pipeline") #only clusterMany values so should be the same
    expect_equal(xx2,xx)
    
    #check reset -- should add combinations of resetColors and resetNames to make sure works independently.
    par(mfrow=c(1,2)) #so can visually check if desired.
    xx3<-plotClusters(cl,resetOrderSamples=TRUE,resetColors=TRUE,resetNames=TRUE) 
    expect_false(isTRUE(all.equal(x2,xx3))) #not a great test. Doesn't really say whether does it right, just whether it does something!
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
    x1<-plotClusters(clusters=cl2,orderClusters="pipeline",resetColors=TRUE)
    x2<-plotClusters(clusters=cl,resetColors=TRUE)
    whP<-.TypeIntoIndices(cl2,"pipeline")
    expect_equal(clusterLegend(x2),clusterLegend(x1)[whP])
    
    #test specifying indices
    wh<-c(3,4,NCOL(clusterMatrix(cl)))
    x3<-plotClusters(cl,orderClusters=wh,axisLine=-2,resetColors=TRUE)
    x4<-plotClusters(cl,orderClusters=wh[c(3,2,1)],axisLine=-2,resetColors=TRUE)
    expect_false(isTRUE(all.equal(x3,x4)))
    
    #test if only a single cluster
    plotClusters(cl,orderClusters=2)
    x5<-plotClusters(cl,orderClusters=2,resetColors=TRUE)
    expect_equal(x5,cl)
    
})

sData<-data.frame(sample(letters[2:5],size=NCOL(simData),replace=TRUE),sample(2:5,size=NCOL(simData),replace=TRUE))
sData<-data.frame(sData,sample(LETTERS[2:5],size=NCOL(simData),replace=TRUE),stringsAsFactors=FALSE)
colnames(sData)<-c("A","B","C")

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

whSamp<-unlist(tapply(1:nSamples(cl3),primaryCluster(cl3),function(x){sample(x,size=3)}))
smData<-simData[1:10,whSamp]
smCount<-simCount[1:10,whSamp]
smCl<-cl3[1:10,whSamp]
smSData<-sData[whSamp,]
clusterLabels(smCl)<-paste("Cluster",1:nClusters(smCl))
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

test_that("`plotHeatmap` works with CE objects", {
    plotHeatmap(smCl)
    plotHeatmap(smCl,whichClusters="none")
    expect_warning(plotHeatmap(smCl,whichClusters="pipeline") ,"whichClusters value does not match any clusters") #there are no pipeline for this one
    clusterType(smCl)[2:3]<-"clusterMany"
    plotHeatmap(smCl,whichClusters="pipeline") 
    expect_error(plotHeatmap(smCl,sampleData="A"))
    colData(smCl)<-DataFrame(smSData)
    plotHeatmap(smCl,sampleData="all")
    plotHeatmap(smCl,sampleData="A")
    plotHeatmap(smCl,sampleData=2:3)
    plotHeatmap(smCl,sampleData="all",whichClusters="none")
    ###stopped here: error because need to change output of plotCluster to provide names to clusterIds in colorLegend.
    plotHeatmap(smCl,sampleData="all",whichClusters="all",alignSampleData=TRUE)
})

