############
## Checks that make up validity of ClusterExperiment:
############


.checkMerge<-function(object){
  ############
  ##Check merge related slots
  ############
	
  if(!is.null(object@dendro_clusters)){
		#further checks with dendro
  	data.cl<-phylobase::tdata(object@dendro_clusters,type="all")
  	if(any(!gsub("ClusterId","",na.omit(data.cl$ClusterIdMerge)) %in% as.character(object@clusterMatrix[,object@merge_index]))) return("ClusterIdMerge information in dendrogram slot must match the corresponding cluster ids in clustering defined by merge_index slot")
	}
	
  #these slots must all be same (either all empty or all set)
  emptyMergeMethod<-is.na(object@merge_method) || length(object@merge_method)==0
  emptyMergeIndex<- is.na(object@merge_index) || length(object@merge_index)==0
  if(emptyMergeIndex & !emptyMergeMethod) return("merge_index NA but merge_method has value")
  if(!emptyMergeIndex & emptyMergeMethod) return("if merge_index not NA, must have value for merge_method")
  emptyCutoff<-is.na(object@merge_cutoff) || length(object@merge_cutoff)==0
  if(emptyMergeIndex & !emptyCutoff) return("merge_index NA but merge_cutoff has value")
  if(!emptyMergeIndex & emptyCutoff) return("if merge_index not NA, must have value for merge_cutoff")
  if(emptyMergeIndex & !is.null(object@merge_nodeMerge)) return("merge_index NA but merge_nodeMerge has value")
  if(!emptyMergeIndex & is.null(object@merge_nodeMerge)) return("if merge_index not NA, must have value for merge_nodeMerge")
  
  #these slots can be set even if no merge, but not vice versa
  emptyDEMethod<-is.na(object@merge_demethod) || length(object@merge_demethod)==0
  if(!emptyMergeIndex & emptyDEMethod) return("if merge_index not NA, must have value for merge_demethod") 
  emptyDendroIndex<-is.na(object@merge_dendrocluster_index) || length(object@merge_dendrocluster_index)==0
  if(!emptyMergeIndex & emptyDendroIndex) return("if merge_index not NA, must have value for merge_dendrocluster_index")
  
  ##Check when there is actual merging:
  if(!emptyMergeIndex){
    if(object@merge_cutoff>1 || object@merge_cutoff<0) return("merge_cutoff should be between 0 and 1")
    if(object@merge_index==object@merge_dendrocluster_index) return("merge_index should not be same as merge_dendrocluster_index")
    if(!length(object@merge_method)==1) return("merge_method must be of length 1")
    
    #deal with possible FC
    baseMergeMethod<-sapply(strsplit(object@merge_method,"_"),.subset2,1)
    if(!baseMergeMethod %in% .availMergeMethods) return(paste("merge_method must be one of available merge methods:", paste(.availMergeMethods,collapse=",")," (with possibility of fold-change added to method name for 'adjP')"))
    allowMergeColumns<-c('Contrast','isMerged','mergeClusterId','NodeId')
    if(!identical(sort(colnames(object@merge_nodeMerge)),sort(allowMergeColumns)) ) {
      return(paste("merge_nodeMerge must have 5 columns and column names equal to:",paste(allowMergeColumns,collapse=",")))
    }
    if(!is.character(object@merge_nodeMerge[,"NodeId"])) return("'NodeId' column of merge_nodeMerge must be character")
    if(!is.character(object@merge_nodeMerge[,"Contrast"])) return("'Contrast' column of merge_nodeMerge must be character")
    if(!is.logical(object@merge_nodeMerge[,"isMerged"])) return("'isMerged' column of merge_nodeMerge must be logical")
    if(!is.numeric(object@merge_nodeMerge[,"mergeClusterId"]) & !all(is.na(object@merge_nodeMerge[,"mergeClusterId"]))) return("'mergeClusterId' column of merge_nodeMerge must be numeric")
    if(any(object@merge_nodeMerge[,"isMerged"])){
      wh<-which(object@merge_nodeMerge[,"isMerged"])
      if(all(is.na(object@merge_nodeMerge[wh,"mergeClusterId"]))) return("mergeClusterId entries of merge_nodeMerge cannot be all NA if there isMerged column that is TRUE")
    }    
    id<-object@merge_nodeMerge[,"mergeClusterId"]
    merg<-object@merge_nodeMerge[,"isMerged"]
    if(length(unique(na.omit(id))) != length(na.omit(id))) return("'mergeClusterId values in merge_nodeMerge not unique")
    if(any(!is.na(id) & !merg)) return("Cannot have values 'mergeClusterId' where 'isMerged' is FALSE")
    cl<-clusterMatrix(object)[,object@merge_index]
    if(any(!na.omit(id)%in% cl)) return("Values in 'mergeClusterId' not match cluster id values")
	
  }
  if(!is.null(object@merge_nodeProp)){
    if(emptyDendroIndex){return("merge_nodeProp is not NULL but merge_dendrocluster_index has no value")
      
    }
    if(length(object@merge_demethod)!=1 || !object@merge_demethod %in% .demethods)
      return(paste("merge_demethod must be one of:",paste(.demethods,collapse=",")))
		
    requireColumns<-c("NodeId","Contrast",.availMergeMethods)
    #need to allow for log fold change columns of adjP
    allCnames<-colnames(object@merge_nodeProp)
    whFC<-grep("adjP_",allCnames)
    whNode<-which(allCnames %in% c("NodeId","Contrast"))
    namesToCheck<-if(length(whFC)>0) allCnames[-whFC] else allCnames
    if(!identical(sort(namesToCheck),sort(requireColumns)) ) 
      return(paste("merge_nodeProp must be data.frame with at least",length(requireColumns),"columns that have column names equal to:",paste(requireColumns,sep="",collapse=",")))
    
    if(!is.character(object@merge_nodeProp[,"NodeId"])) return("'Node' column of merge_nodeProp must be character")
    if(!is.character(object@merge_nodeProp[,"Contrast"])) return("'Contrast' column of merge_nodeProp must be character")
    for(method in allCnames[-whNode]){
      if(!is.numeric(object@merge_nodeProp[,method])) return(paste(method,"column of merge_nodeProp must be numeric"))
    }
	
  }
  #Check that dendro node ids match those in merge tables.
  emptyDendro<-is.na(object@dendro_index) || length(object@dendro_index)==0
  if(!emptyDendroIndex && !emptyDendro && object@merge_dendrocluster_index==object@dendro_index){
	dendroNodes<-phylobase::tdata(object@dendro_clusters,type="internal")[,"NodeId"]
	nodes<-object@merge_nodeMerge[,"NodeId"]
	if(!all(sort(nodes) == sort(dendroNodes))) return("Not all of nodes in dendro_clusters have a value in merge_nodeMerge")
	nodes<-object@merge_nodeProp[,"NodeId"]
	if(!all(sort(nodes) == sort(dendroNodes))) return("Not all of nodes in dendro_clusters have a value in merge_nodeProp")		
  }
  return(TRUE)
}
.checkCoClustering<-function(object){
    ## FIXME: this needs to be changed if now return NxB matrix. 
    ## For now, commenting out the check entirely. Can be anything...
    ## Check co-clustering
    # if(!is.null(object@coClustering) &&
    #    (NROW(object@coClustering) != NCOL(object@coClustering)
    #     | NCOL(object@coClustering) != NCOL(object))) {
    #     return("`coClustering` must be a sample by sample matrix.")
    # }
    typeCoCl<-.typeOfCoClustering(object)
    if(typeCoCl=="indices"){
        if(!all(object@coClustering %in% 1:nClusterings(object)) ) return("CoClustering slot is a vector, but doesn't match indices of clusterMatrix of the object")
    }
    return(TRUE)
}


#' @rdname ClusterFunction-class
#' @export
#' @aliases internalFunctionCheck
#' @examples
#' #Use internalFunctionCheck to check possible function
#' goodFUN<-function(inputMatrix,k,cluster.only,...){
#'	cluster::pam(x=t(inputMatrix),k=k,cluster.only=cluster.only)
#' }
#' #passes internal check
#' internalFunctionCheck(goodFUN,inputType=c("X","diss"),
#'    algorithmType="K",outputType="vector")
#' myCF<-ClusterFunction(clusterFUN=goodFUN, inputType="X",
#'    algorithmType="K", outputType="vector")
#' #doesn't work, because haven't made results return vector when cluster.only=TRUE
#' badFUN<-function(inputMatrix,k,cluster.only,...){cluster::pam(x=inputMatrix,k=k)}
#' internalFunctionCheck(badFUN,inputType=c("X","diss"),
#'    algorithmType="K",outputType="vector")
#' @details \code{internalFunctionCheck} is the function that is called by the 
#'   validity check of the \code{ClusterFunction} constructor (if 
#'   \code{checkFunctions=TRUE}). It is available as an S3 function for the user
#'   to be able to test their functions and debug them, which is difficult to do
#'   with a S4 validity function.
#' @return Returns a logical value of TRUE if there are no problems. If there is 
#'   a problem, returns a character string describing the problem encountered.
internalFunctionCheck<-function(clusterFUN,inputType,algorithmType,outputType){
    #--- Make small data
    N<-20
    #--- Set parameters
    if(algorithmType=="01") argList<-list(alpha=.5)	
    if(algorithmType=="K") argList<-list(k=2)	
    argList<-c(argList,list(cluster.only=TRUE,checkArgs=FALSE))
    #--- Run function on small data
    if("X" %in% inputType){
        set.seed(2851)
        x<-matrix(rnorm(N*3),ncol=N,nrow=3)
        test<-try(do.call(clusterFUN, 
            c(list(inputMatrix=x,inputType="X"), argList)),silent=TRUE)
        if(inherits(test,"try-error")) return(paste("function test fails with inputType='X'. ",test[1]))
    }
    if("diss" %in% inputType){
        set.seed(2851)
        diss<-matrix(runif(N^2,min=0,max=0.5),ncol=N,nrow=N)
        diss<-diss + t(diss)
        diag(diss)<-0
        test<-try(do.call(clusterFUN,
            c(list(inputMatrix=diss,inputType="diss"),argList)),silent=TRUE)
        if(inherits(test,"try-error")) return(paste("function test fails with input inputType='diss'.",test[1]))
    }
    if("cat" %in% inputType){
        set.seed(2851)
        x<-matrix(sample(x=1:4,size=N*3,replace=TRUE),ncol=N,nrow=3)
        x<-cbind(x,x)
        x<-x[,sample(1:ncol(x))]
        test<-try(do.call(clusterFUN,
            c(list(inputMatrix=x,inputType="cat"),argList)),silent=TRUE)
        if(inherits(test,"try-error")) return(paste("function test fails with inputType='cat'. ",test[1]))
    }
    if(outputType=="vector"){
        if(length(test)!=N) return("clusterFUN does not return a vector equal to the number of observations")
    }
    return(TRUE)
}


.checkHasArgs<-function(FUN,requiredArgs){
    funArgs<-names(as.list(args(FUN)))
    all(requiredArgs %in% funArgs)
}
