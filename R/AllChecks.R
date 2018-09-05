############
## Checks that make up validity of ClusterExperiment:
############
.checkAssays<-function(object){
  # ############
  # check assays
  # ############
  if(length(assays(object)) < 1) {
    return("There must be at least one assay slot.")
  }
  if(!is.numeric(assay(object))) {
    return("The data must be numeric.")
  }
  if(anyNA(assay(object))) {
    return("NA values are not allowed.")
  }
  return(TRUE)
}
.checkTransform<-function(object){
  #############
  #check transform
  #############
  tX <- try(transformData(object),silent=TRUE)
  if(inherits(tX, "try-error")){
    stop(paste("User-supplied `transformation` produces error on the input data
                 matrix:\n",tX))
  }
  if(anyNA(tX)) {
    return("NA values after transforming data matrix are not allowed.")
  }
  return(TRUE)
}

.checkClusterMatrix<-function(object){
  ############
  #Check clusterMatrix
  ############
  zeroRow<-NROW(object@clusterMatrix)==0 #could happen in subsetting
  if(!all(is.na((object@clusterMatrix))) &
     !(NROW(object@clusterMatrix) == NCOL(object))) {
    return("If present, `clusterMatrix` must have as many row as cells.")
  }
  if(!zeroRow & !is.numeric(object@clusterMatrix)) {
    return("`clusterMatrix` must be a numeric matrix.")
  }
  
  if(NCOL(object@clusterMatrix)!= length(object@clusterTypes)) {
    return("length of clusterTypes must be same as NCOL of the clusterMatrix")
  }
  
  if(NCOL(object@clusterMatrix)!= length(object@clusterInfo)) {
    return("length of clusterInfo must be same as NCOL of the clusterMatrix")
  }
  #check internally stored as integers
  testConsecIntFun<-function(x){
    whCl<-which(!x %in% c(-1,-2))
    uniqVals<-unique(x[whCl])
    return(all(sort(unname(uniqVals))==seq_along(uniqVals)))
  }
	if(!zeroRow){
	  testConsecIntegers<-apply(object@clusterMatrix,2,testConsecIntFun)
	  if(!all(testConsecIntegers)) return("the cluster ids in clusterMatrix must be stored internally as consecutive integer values")		
	}
  
  return(TRUE)
}
.checkClusterLabels<-function(object){
  if(!all(is.na(object@clusterMatrix))){
    if(is.null(colnames(object@clusterMatrix))) return("must have clusterLabels by assignment of column names for clusterMatrix")
    if(any(duplicated(colnames(object@clusterMatrix)))) return("cannot have duplicated clusterLabels (i.e. clusterMatrix must have unique column names)")
    
  }
  return(TRUE)
}
.checkDendrogram<-function(object){
  ############
  ##Check dendrogram
  ############
  if(!is.null(object@dendro_samples)){
    if(nobs(object@dendro_samples) != NCOL(object)) {
      return("dendro_samples must have the same number of leaves as the number of samples")
    }
    if(is.na(object@dendro_outbranch)) return("if dendro_samples is defined, must also define dendro_outbranch")
  }
  else{
    # if(!is.null(object@dendro_clusters)) return("dendro_samples should not be null if dendro_clusters is non-null") #now optional to have samples
    if(!is.na(object@dendro_outbranch)) return("dendro_samples should not be null if dendro_outbranch is not NA")
  }
  if(!is.null(object@dendro_clusters)){
    if(is.na(dendroClusterIndex(object))) return("if dendrogram slots are filled, must have corresponding dendro_index defined.")
    dcluster<-clusterMatrix(object)[,dendroClusterIndex(object)]
    if(nobs(object@dendro_clusters) != max(dcluster)) {
      return("dendro_clusters must have the same number of leaves as the number of (non-negative) clusters")
    }
  }
  else{
    # if(!is.null(object@dendro_samples)) return("dendro_clusters should not be null if dendro_samples is non-null") #now optional to have samples
  }
  return(TRUE)
}


.checkMerge<-function(object){
  ############
  ##Check merge related slots
  ############
  #these slots must all be same (either all empty or all set)	
  if(is.na(object@merge_index) & !is.na(object@merge_method)) return("merge_index NA but merge_method has value")
  if(!is.na(object@merge_index) & is.na(object@merge_method)) return("if merge_index not NA, must have value for merge_method")
  if(is.na(object@merge_index) & !is.na(object@merge_cutoff)) return("merge_index NA but merge_cutoff has value")
  if(!is.na(object@merge_index) & is.na(object@merge_cutoff)) return("if merge_index not NA, must have value for merge_cutoff")
  if(is.na(object@merge_index) & !is.null(object@merge_nodeMerge)) return("merge_index NA but merge_nodeMerge has value")
  if(!is.na(object@merge_index) & is.null(object@merge_nodeMerge)) return("if merge_index not NA, must have value for merge_nodeMerge")
  
  #these slots can be set even if no merge, but not vice versa
  if(!is.na(object@merge_index) & is.na(object@merge_demethod)) return("if merge_index not NA, must have value for merge_demethod")  
  if(!is.na(object@merge_index) & is.na(object@merge_dendrocluster_index)) return("if merge_index not NA, must have value for merge_dendrocluster_index")
  
  ##Check when there is actual merging:
  if(!is.na(object@merge_index)){
    if(object@merge_cutoff>1 || object@merge_cutoff<0) return("merge_cutoff should be between 0 and 1")
    if(object@merge_index==object@merge_dendrocluster_index) return("merge_index should not be same as merge_dendrocluster_index")
    if(!length(object@merge_method)==1) return("merge_method must be of length 1")
    
    #deal with possible FC
    baseMergeMethod<-sapply(strsplit(object@merge_method,"_"),.subset2,1)
    if(!baseMergeMethod %in% .availMergeMethods) return(paste("merge_method must be one of available merge methods:", paste(.availMergeMethods,collapse=",")," (with possibility of fold-change added to method name for 'adjP')"))
    allowMergeColumns<-c('Contrast','isMerged','mergeClusterId','Node')
    
    if(!identical(sort(colnames(object@merge_nodeMerge)),sort(allowMergeColumns)) ) {
      return(paste("merge_nodeMerge must have 4 columns and column names equal to:",paste(allowMergeColumns,collapse=",")))
    }
    if(!is.character(object@merge_nodeMerge[,"Node"])) return("'Node' column of merge_nodeMerge must be character")
    if(!is.character(object@merge_nodeMerge[,"Contrast"])) return("'Contrast' column of merge_nodeMerge must be character")
    if(!is.logical(object@merge_nodeMerge[,"isMerged"])) return("'isMerged' column of merge_nodeMerge must be logical")
    if(!is.numeric(object@merge_nodeMerge[,"mergeClusterId"]) & !all(is.na(object@merge_nodeMerge[,"mergeClusterId"]))) return("'mergeClusterId' column of merge_nodeMerge must be numeric")
    if(any(object@merge_nodeMerge[,"isMerged"])){
      wh<-which(object@merge_nodeMerge[,"isMerged"])
      if(all(is.na(object@merge_nodeMerge[wh,"mergeClusterId"]))) return("mergeClusterId entries of merge_nodeMerge cannot be all NA if there isMerged column that is TRUE")
    }    
    id<-    object@merge_nodeMerge[,"mergeClusterId"]
    merg<-object@merge_nodeMerge[,"isMerged"]
    if(length(unique(na.omit(id))) != length(na.omit(id))) return("'mergeClusterId values in merge_nodeMerge not unique")
    if(any(!is.na(id) & !merg)) return("Cannot have values 'mergeClusterId' where 'isMerged' is FALSE")
    cl<-clusterMatrix(object)[,object@merge_index]
    if(any(!na.omit(id)%in% cl)) return("Values in 'mergeClusterId' not match cluster id values")
  }
  if(!is.null(object@merge_nodeProp)){
    if(is.na(object@merge_dendrocluster_index)){return("merge_nodeProp is NULL but merge_dendrocluster_index has value")
      
    }
    if(length(object@merge_demethod)!=1 || !object@merge_demethod %in% .demethods)
      return(paste("merge_demethod must be one of:",paste(.demethods,collapse=",")))
		
    requireColumns<-c("Node","Contrast",.availMergeMethods)
    #need to allow for log fold change columns of adjP
    allCnames<-colnames(object@merge_nodeProp)
    whFC<-grep("adjP_",allCnames)
    whNode<-which(allCnames %in% c("Node","Contrast"))
    namesToCheck<-if(length(whFC)>0) allCnames[-whFC] else allCnames
    if(!identical(sort(namesToCheck),sort(requireColumns)) ) 
      return(paste("merge_nodeProp must be data.frame with at least",length(requireColumns),"columns that have column names equal to:",paste(requireColumns,sep="",collapse=",")))
    
    if(!is.character(object@merge_nodeProp[,"Node"])) return("'Node' column of merge_nodeProp must be character")
    if(!is.character(object@merge_nodeProp[,"Contrast"])) return("'Contrast' column of merge_nodeProp must be character")
    for(method in allCnames[-whNode]){
      if(!is.numeric(object@merge_nodeProp[,method])) return(paste(method,"column of merge_nodeProp must be numeric"))
    }
  }
  return(TRUE)
}
.checkCoClustering<-function(object){
  ## Check co-clustering
  if(!is.null(object@coClustering) &&
     (NROW(object@coClustering) != NCOL(object@coClustering)
      | NCOL(object@coClustering) != NCOL(object))) {
    return("`coClustering` must be a sample by sample matrix.")
  }
  return(TRUE)
}

.checkPrimaryIndex<-function(object){
  ############
  ## Check primary index
  ############
  if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
    #check primary index
    if(length(object@primaryIndex) != 1) {
      if(length(object@primaryIndex) == 0) return("If more than one set of clusterings, a primary cluster must be specified.")
      if(length(object@primaryIndex) > 0) return("Only a single primary index may be specified")
    }
    if(object@primaryIndex > NCOL(object@clusterMatrix) | object@primaryIndex < 1) {
      return("`primaryIndex` out of bounds.")
    }
  }
  return(TRUE)
}

.checkClusterTypes<-function(object){
  ############
  ## Check clusterTypes
  ############
  if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
    if(NCOL(object@clusterMatrix) != length(object@clusterTypes)) {
      return("`clusterTypes` must be the same length as NCOL of `clusterMatrix`.")
    }
    if(!is.null(names(object@clusterTypes))) return("clusterTypes should not have names")
  }
  return(TRUE)
}

#this check is just to check object, not whether matches clusters
.checkClusterLegendList<-function(clusterLegend,allowNames=TRUE,reqNames=c("clusterIds","color","name")){
  if(!is.list(clusterLegend)) return("clusterLegend must be a list")
  if(!allowNames & !is.null(names(clusterLegend))) return("clusterLegend should not have names")
  testIsMatrix <- sapply(clusterLegend, function(x) {!is.null(dim(x))})
  if(!all(testIsMatrix)) {
    return("Each element of `clusterLegend` list must be a matrix")
  }
  testNames<-sapply(clusterLegend,function(x){
    if(is.null(colnames(x))) return(FALSE)
    else{
      if(!all(reqNames %in% colnames(x))) return(FALSE)
      else return(TRUE)
    }
  })
  if(!all(testNames)) {
    return(paste("each element of `clusterLegend` must be matrix with column names defined, with at a minimum the names names:", paste(reqNames,collapse=",")))
  }
  testColorCols1 <- sapply(clusterLegend, function(x){is.character(x)})
  if(!all(testColorCols1)) {
    return("each element of `clusterLegend` must be matrix of character values")
  }
  return(TRUE)
}
#this check checks both object (calls .checkClusterLegend) and whether matches clusters
#make so can call on arbitrary clusterLegend...not need to be CE object
.checkClustersWithClusterLegend<-function(clusters,clusterLegend){
  #check structure of clusterLegend list -- for CE object, can't have names.
  legendCheck<-.checkClusterLegendList(clusterLegend,allowNames=FALSE,reqNames=c("clusterIds","color","name"))
  if(!is.logical(legendCheck)) return(legendCheck)
  
  #check matches clusters
  if(length(clusterLegend) != NCOL(clusters)) {
    return("`clusterLegend` must be list of same length as NCOL of
               `clusterMatrix`")
  }
  testColorRows <- sapply(clusterLegend, function(x){nrow(x)})
  testClusterMat <- apply(clusters, 2, function(x) {length(unique(x))})
  if(!all(testColorRows == testClusterMat)) {
    return("each element of `clusterLegend` must be matrix with number of rows equal to the number of clusters (including -1 or -2 values) in `clusterMatrix`")
  }
  testColorCols1 <- sapply(seq_along(clusterLegend), function(ii){
    col<-clusterLegend[[ii]]
    x<-clusters[,ii]
    y<-col[,"clusterIds"]
    if(is.numeric(x)){
      y<-as.numeric(y)
    }
    return(all(y %in% x) & all(x %in% y))
  })
  if( !all(testColorCols1)) {
    return("each element of `clusterLegend` must be matrix with column
             `clusterIds` matching the corresponding
             clusterMatrix values")
  }
  
  return(TRUE)
}
.checkClusterLegend<-function(object){
  if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
    return(.checkClustersWithClusterLegend(clusters=object@clusterMatrix,clusterLegend=object@clusterLegend))
  }
  return(TRUE)
}

.checkOrderSamples<-function(object){
  ####
  #test orderSamples
  ####
  if(length(object@orderSamples)!=NCOL(assay(object))) {
    return("`orderSamples` must be of same length as number of samples
	       (NCOL(assay(object)))")
  }
  if(any(!object@orderSamples %in% seq_len(NCOL(assay(object))))) {
    return("`orderSamples` must be values between 1 and the number of samples.")
  }
  return(TRUE)
}

#' @rdname ClusterFunction-class
#' @export
#' @aliases internalFunctionCheck
#' @examples
#' #Use internalFunctionCheck to check possible function
#' goodFUN<-function(x,diss,k,checkArgs,cluster.only,...){
#'	cluster::pam(x=t(x),k=k,cluster.only=cluster.only)
#' }
#' #passes internal check
#' internalFunctionCheck(goodFUN,inputType="X",algorithmType="K",outputType="vector")
#' #Note it doesn't pass if inputType="either" because no catches for x=NULL
#' internalFunctionCheck(goodFUN, inputType="either",algorithmType="K",outputType="vector")
#' myCF<-ClusterFunction(clusterFUN=goodFUN, inputType="X",algorithmType="K", outputType="vector")
#' badFUN<-function(x,diss,k,checkArgs,cluster.only,...){cluster::pam(x=x,k=k)}
#' internalFunctionCheck(badFUN,inputType="X",algorithmType="K",outputType="vector")
#' @details \code{internalFunctionCheck} is the function that is called by the 
#'   validity check of the \code{ClusterFunction} constructor (if 
#'   \code{checkFunctions=TRUE}). It is available as an S3 function for the user
#'   to be able to test their functions and debug them, which is difficult to do
#'   with a S4 validity function.
internalFunctionCheck<-function(clusterFUN,inputType,algorithmType,outputType){
  #--- Make small data
  N<-20
  set.seed(2851)
  x<-matrix(rnorm(N*3),ncol=N,nrow=3)
  set.seed(2851)
  diss<-matrix(runif(N^2,min=0,max=0.5),ncol=N,nrow=N)
  diss<-diss + t(diss)
  diag(diss)<-0
  #--- Set parameters
  if(algorithmType=="01") argList<-list(alpha=.5)	
  if(algorithmType=="K") argList<-list(k=2)	
  argList<-c(argList,list(cluster.only=TRUE,checkArgs=FALSE))
  #--- Run function on small data
  if(inputType %in% c("X")){
    test<-try(do.call(clusterFUN,c(list(x=x),argList)),silent=TRUE)
    if(inherits(test,"try-error")) return(paste("function test fails with input X. ",test[1]))
  }
  if(inputType %in% c("diss")){
    test<-try(do.call(clusterFUN,c(list(diss=diss),argList)),silent=TRUE)
    if(inherits(test,"try-error")) return(paste("function test fails with input diss.",test[1]))
  }
  if(inputType %in% c("either")){
    test1<-try(do.call(clusterFUN,c(list(x=x,diss=NULL),argList)),silent=TRUE)
    if(inherits(test1,"try-error")) return(paste("function test fails with input x and NULL diss.",test1[1]))
    test2<-try(do.call(clusterFUN,c(list(x=NULL,diss=diss),argList)),silent=TRUE)
    if(inherits(test2,"try-error")){
      return(paste("function test fails with input diss and NULL x.",test2[1]))
    }
    test3<-try(do.call(clusterFUN,c(list(x=x,diss=diss),argList)),silent=TRUE)
    if(inherits(test3,"try-error")) return(paste("function test fails both diss and x input.",test3[1]))
    if(outputType=="vector" & length(test1)!=N || length(test2)!=N || length(test3)!=N) return("clusterFUN does not return a vector equal to the number of observations.")
  }
  else{
    if(outputType=="vector"){
      if(length(test)!=N) return("clusterFUN does not return a vector equal to the number of observations")
    }
  }
  return(TRUE)
}


.checkHasArgs<-function(FUN,requiredArgs){
  funArgs<-names(as.list(args(FUN)))
  all(requiredArgs %in% funArgs)
}