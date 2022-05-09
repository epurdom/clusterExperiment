.mynote<-function(x){
	ClusterExperiment:::.mynote
}
# Small function to identify what type of CoClustering information is stored in co-clustering slot
# (might make sense to make it a slot, but then have to update object, for something not so important...)
.typeOfCoClustering<-function(ceObj){
	if(is.null(ceObj@coClustering)) return("null")
    if(is.null(dim(ceObj@coClustering))) return("indices")
        else if(!isSymmetric(ceObj@coClustering)) return("clusterings")
            else return("distance")
}

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
    retval<-clusterExperiment:::.addBackSEInfo(newObj,oldObj) #Creates ClusterExperiment object with everything. 
    retval<-RSECClass(retval,
                merge_index=newObj@merge_index,
                merge_cutoff=newObj@merge_cutoff,
                merge_dendrocluster_index=newObj@merge_dendrocluster_index,
                merge_nodeProp=newObj@merge_nodeProp,
                merge_nodeMerge=newObj@merge_nodeMerge,
                merge_method=newObj@merge_method,
                merge_demethod=newObj@merge_demethod,
                coClustering=coClustering(newObj)
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