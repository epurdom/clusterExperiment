
makeBlankData <- function(data,groupsOfFeatures,nBlankLines = 1) {
  if (!is.list(groupsOfFeatures))
    stop("groupsOfFeatures must be a list. Will ignore.")
  
  testIndices <- sapply(groupsOfFeatures,function(x) {
    is.numeric(x) && all(x %in% 1:NROW(data))
  })
  if (!all(testIndices))
    stop("Invalid list of indices in groupsOfFeatures.")
  if(is.null(rownames(data))) row.names(data)<-as.character(1:nrow(data))
  #make list of data of feature groups
  dataList <- lapply(groupsOfFeatures,function(ii){data[ii,,drop=FALSE]})

  #add NA rows between groups
  naData <- matrix(NA,nrow = nBlankLines,ncol = ncol(data))
  dataListMinus <- lapply(dataList[-length(dataList)],function(x) {
    return(rbind(x,naData))
  })
  newData <-
    data.frame(do.call("rbind",c(dataListMinus,dataList[length(dataList)])))
  
  #make names for this
  rnames <- lapply(dataList,rownames)
  rnamesMinus <-
    lapply(head(rnames,-1),function(x) {
      c(x,rep("",nBlankLines))
    })
  rnames <- unname(c(unlist(rnamesMinus),rnames[[length(rnames)]]))
  #browser()
  #can't set rownames because not unique values!
  #rownames(newData)<-rnames
  return(list(dataWBlanks = newData,rowNamesWBlanks = rnames))
  
}
# temp<-makeBlankData(norm[geneord,keep],sep=head(sep,-1),nBlankLines=2)
# dat<-t(temp$data)
# rnames<-temp$rownames
# 
# out <- dualHeatmap(clusterVec=cl, heatData=t(norm[,keep]), dual=FALSE, clusterSamples=FALSE, clusterVars=FALSE, annCol=data.frame(Cluster=cl), annColors=list(Cluster=colMerged[c(7, 2, 5, 8)]), main="", whVars = geneord, annLegend=FALSE)
# 
# dualHeatmap(clusterVec=cl, heatData=dat, dual=FALSE, clusterSamples=FALSE, clusterVars=FALSE, annCol=data.frame(Cluster=cl), annColors=list(Cluster=colMerged[c(7, 2, 5, 8, 6)]), main="Allen Institute L5 markers", labRow=rnames, annLegend=FALSE, breaks=out$breaks)
