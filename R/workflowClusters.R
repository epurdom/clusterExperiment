#Update here if change workflow values. Also defines the order of them.
.workflowValues<-c("final","mergeClusters","combineMany","clusterMany")

#' @title Methods for workflow clusters
#'
#' @description The main workflow of the package is made of
#'   \code{\link{clusterMany}}, \code{\link{combineMany}}, and
#'   \code{\link{mergeClusters}}. The clusterings from these functions (and not
#'   those obtained in a different way) can be obtained with the functions
#'   documented here.
#' @param x a \code{\link{ClusterExperiment}} object.
#' @param iteration numeric. Which iteration of the workflow should be used.
#' @return \code{workflowClusters} returns a matrix consisting of the
#'   appropriate columns of the \code{clusterMatrix} slot.
#' @aliases workflowClusters workflowClusterTable workflowClusterDetails
#'   workflowClusters,ClusterExperiment-method
#' @name workflowClusters
#' @export
#' @examples
#' data(simData)
#'
#' cl <- clusterMany(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(FALSE), removeSil=TRUE,
#' subsample=FALSE)
#'
#' clCommon <- combineMany(cl, whichClusters="workflow", proportion=0.7,
#' minSize=10)
#'
#' clCommon <- makeDendrogram(clCommon)
#'
#' clMerged <- mergeClusters(clCommon,mergeMethod="adjP")
#'
#' head(workflowClusters(clMerged))
#' workflowClusterDetails(clMerged)
#' workflowClusterTable(clMerged)
setMethod(
  f = "workflowClusters",
  signature = signature("ClusterExperiment"),
  definition = function(x,iteration=0) {
    ppIndex<-workflowClusterDetails(x)
    if(is.na(iteration)) iteration<-unique(ppIndex[,"iteration"])
    if(!is.null(ppIndex)){
      whIteration<-which(ppIndex[,"iteration"]%in%iteration)
      if(length(whIteration)>0){
        index<-ppIndex[whIteration,"index"]
        return(clusterMatrix(x)[,index,drop=FALSE])
      }
      else return(NULL)
    }
    else return(NULL)
  }
)


#' @rdname workflowClusters
#' @return \code{workflowClusterDetails} returns a \code{data.frame} with some
#'   details on the clusterings, such as the type (e.g., `clusterMany`,
#'   `combineMany`) and iteration.
#' @export
setMethod(
  f = "workflowClusterDetails",
  signature = signature("ClusterExperiment"),
  definition = function(x) {

    if(length(clusterType(x))!=NCOL(clusterMatrix(x))) stop("Invalid ClusterExperiment object")
    #check if old iterations already exist; note assumes won't have previous iteration unless have current one.
    existingOld<-lapply(.workflowValues,function(ch){
      regex<-paste(ch,"_",sep="")
      grep(regex,clusterType(x))

    })
    st<-strsplit(clusterType(x)[unlist(existingOld)],"_")
    oldValues<-data.frame(index=unlist(existingOld),type=sapply(st,.subset2,1),iteration=as.numeric(sapply(st,.subset2,2)),stringsAsFactors=FALSE)

    wh<-which(clusterType(x) %in% .workflowValues) #current iteration
    if(length(wh)>0){
      existingValues<-data.frame(index=wh,type=clusterType(x)[wh], iteration=0,stringsAsFactors=FALSE) #0 indicates current iteration
      if(nrow(oldValues)>0) existingValues<-rbind(oldValues,existingValues)
    }
    else{
      if(nrow(oldValues)>0) existingValues<-oldValues
      else   return(NULL)
    }

    return(existingValues)

  }
)
#' @rdname workflowClusters
#' @return \code{workflowClusterTable} returns a table of how many of the
#'   clusterings belong to each of the following possible values: `final`,
#'   `mergeClusters`, `combineMany` and `clusterMany`.
#' @export
setMethod(
  f = "workflowClusterTable",
  signature = signature("ClusterExperiment"),
  definition = function(x){
    ppIndex<-workflowClusterDetails(x)
    table(Type=factor(ppIndex[,"type"],levels=.workflowValues),Iteration=factor(ppIndex[,"iteration"]))
  }
)

#change current workflow to old iteration 
# add number to it if eraseOld=FALSE
# delete ALL workflow if eraseOld=TRUE (not just the current iteration)
.updateCurrentWorkflow<-function(x,eraseOld,newToAdd){
    ppIndex<-workflowClusterDetails(x)
    if(!newToAdd %in% .workflowValues) stop("error in internal coding: newToAdd must be one of .workflowValues. Contact mantainer.")
    whNew<-match(newToAdd, .workflowValues)
    downstreamType<-.workflowValues[1:whNew]
    newX<-x
    if(!is.null(ppIndex)){ #there are pre-existing workflow results
        curr<-ppIndex[ppIndex[,"iteration"]==0,]
        if(any(curr[,"type"] %in% downstreamType) || any(ppIndex[,"iteration"]!=0)){
            if(eraseOld){ 
                #removes all past iterations, not just current, except for current iteration that upstream of new one
                whRm<- union(curr[curr[,"type"] %in% downstreamType, "index"],ppIndex[ppIndex[,"iteration"]!=0,"index"])
                if(length(whRm)==nClusters(x)) return(NULL)
                else newX<-removeClusters(x,whRm) 
            }
            else{
                #otherwise, only current downstream ones exist
                if(any(curr[,"type"] %in% downstreamType)){
                    newIteration<-max(ppIndex[,"iteration"])+1
                    whFix<-curr[curr[,"type"] %in% downstreamType, "index"]
                    updateCluster<-clusterType(x)
                    updateCluster[whFix]<-paste(updateCluster[whFix],newIteration,sep="_")
                    newX@clusterType<-updateCluster          
                }
            }
        }
    }
    newX<-.unnameClusterSlots(newX)
    validObject(newX)
    return(newX)
}
