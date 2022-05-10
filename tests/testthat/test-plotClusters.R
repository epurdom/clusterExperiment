context("plotClusters")



test_that("`plotClusters` works with matrix, ClusterExperiment objects", {

    #test matrix version
    expect_silent(x<-plotClusters(object=clusterMatrix(ceSimCount)))
    expect_equal(dim(clusterMatrix(ceSimCount)),dim(x$colors))
    expect_equal(dim(clusterMatrix(ceSimCount)),dim(x$aligned))
    expect_equal(length(x$clusterLegend),
      ncol(clusterMatrix(ceSimCount)))

    #test CE version
    expect_silent(x<-plotClusters(ceSimCount))
    expect_is(x,"ClusterExperiment")
    expect_equal( x,ceSimCount)


    #check reset -- should add combinations of resetColors and resetNames to make sure works independently.
    par(mfrow=c(1,2)) #so can visually check if desired.
    expect_silent(xx3<-plotClusters(ceSimCount,
      resetOrderSamples=TRUE,resetColors=TRUE,resetNames=TRUE))
    expect_silent(plotClusters(xx3,existingColors="all"))


    nm<-as.numeric(unlist(lapply(clusterLegend(xx3),
      function(x){x[,"name"]})))
    col<-(unlist(lapply(clusterLegend(xx3),
      function(x){x[,"color"]})))
    wh<-which(col %in% c("white","grey"))
    expect_equal(match(col[-wh],bigPalette),nm[-wh])
    nmOld<-as.numeric(unlist(lapply(clusterLegend(ceSimCount),
      function(x){x[,"name"]})))
    expect_false(isTRUE(all.equal(nm,nmOld)))
    idOld<-as.numeric(unlist(lapply(clusterLegend(ceSimCount),
      function(x){x[,"clusterIds"]})))
    idNew<-as.numeric(unlist(lapply(clusterLegend(xx3),
      function(x){x[,"clusterIds"]})))
    expect_equal(idOld,idNew)

    #check existing colors
    x2<-plotClusters(ceSimCount,existingColors="all")

    #test -1
    plotClusters(ceSimCount)


    #test specifying indices
    wh<-c(3,4,NCOL(clusterMatrix(ceSimCount)))
    x3<-plotClusters(ceSimCount,whichClusters=wh,
      axisLine=-2,resetColors=TRUE)
    x4<-plotClusters(ceSimCount,whichClusters=wh[c(3,2,1)],
      axisLine=-2,resetColors=TRUE)
    expect_false(isTRUE(all.equal(x3,x4)))

    par(mfrow=c(1,1)) #otherwise will affect other tests.
})



test_that("`plotClusters` rerun above tests with colData included", {

    #test matrix version
    expect_silent(x<-plotClusters(object=clusterMatrix(ceSimCount),
        colData=as.data.frame(colData(ceSimCount))))
    expect_equal(ncol(clusterMatrix(ceSimCount))+
      ncol(colData(ceSimCount)), ncol(x$colors))
    expect_equal(ncol(clusterMatrix(ceSimCount))+
      ncol(colData(ceSimCount)),ncol(x$aligned))
    expect_equal(length(x$clusterLegend),
        ncol(clusterMatrix(ceSimCount))+ncol(colData(ceSimCount)))

    #---    
    #test CE version
    #---    
    
    expect_error(plotClusters(ceSimCount,
        colData=as.data.frame(colData(ceSimCount))),
        "invalid values for pulling sample data from colData of object")
    expect_silent(plotClusters(ceSimCount,colData="all"))
    par(mfrow=c(1,2))
    expect_silent(x2<-plotClusters(ceSimCount,colData="all",resetColors=TRUE))
    expect_silent(x1<-plotClusters(ceSimCount,resetColors=TRUE))


    #check NAs
    naSim<-ceSimCount
    colData(naSim)[sample(size=10,x=1:nrow(naSim)),c("A","B")]<-NA
    expect_silent(plotClusters(naSim,colData=c("A","B")))

    #test the new TRUE option for colData
    expect_silent(plotClusters(naSim,colData=TRUE))

  #this is not working because first one puts -1/-2 last and second puts them first, and so then assigns different colors to the groups
#  expect_equal(x1,x2)
#   par(mfrow=c(1,2))
#   x2<-plotClusters(ceSimCount,colData="all",resetColors=FALSE)
#   x1<-plotClusters(ceSimCount,resetColors=FALSE)
    par(mfrow=c(1,1))

})


