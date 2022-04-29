.mynote<-function(x){
	ClusterExperiment:::.mynote
}
.typeOfCoClustering<-function(x){ClusterExperiment:::.typeOfCoClustering}
## FIXME: Some day, if update addClusterings so could choose where to put new clusterings, could just use that so don't replicate code here.
.addNewResult<-function(newObj,oldObj){
    # want most recent addition on top of clusterMatrix
    # mergeObjects=TRUE means will keep any CE info in newObj, but add CE info in oldObj if not in new object (used to be done by this function, but updated `addClusterings` to do that for me so all in one place). 
    retval<-addClusterings(x=newObj, y= oldObj, 
        transferFrom="x", mergeObjects=TRUE) 
    retval<-.addBackSEInfo(newObj=retval,oldObj=oldObj) #make sure keeps SE info
    #   Note: .addBackSEInfo calls ClusterExperiment (i.e. validates)
    return(retval)
}

#this function keeps everything from new, except grabs SE info from old
.addBackSEInfo<-function(newObj,oldObj){
  retval<-ClusterExperiment(as(oldObj,"SingleCellExperiment"),
                            clusters=clusterMatrix(newObj),
                            transformation=transformation(newObj),
                            clusterTypes=clusterTypes(newObj),
                            clusterInfo=clusteringInfo(newObj),
                            orderSamples=orderSamples(newObj),
                            coClustering=coClustering(newObj),
                            dendro_samples=newObj@dendro_samples,
                            dendro_clusters=newObj@dendro_clusters,
                            dendro_index=newObj@dendro_index,
                            merge_index=newObj@merge_index,
                            merge_cutoff=newObj@merge_cutoff,
                            merge_dendrocluster_index=newObj@merge_dendrocluster_index,
                            merge_nodeProp=newObj@merge_nodeProp,
                            merge_nodeMerge=newObj@merge_nodeMerge,
                            merge_method=newObj@merge_method,
                            merge_demethod=newObj@merge_demethod,
                            primaryIndex=primaryClusterIndex(newObj),
                            clusterLegend=clusterLegend(newObj),
                            checkTransformAndAssay=FALSE
  )
  return(retval)
}
.addPrefixToClusterNames<-function(ceObj,prefix,whCluster){
  ceLegend<-clusterLegend(ceObj)[[whCluster]]
  cl<-ceLegend[,"clusterIds"]
  ceLegend[,"name"]<-ClusterExperiment::numericalAsCharacter(values=cl,prefix=prefix)
  clusterLegend(ceObj)[[whCluster]]<-ceLegend
  return(ceObj)
}