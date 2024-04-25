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

#these functions are checks where don't need the corresponding object information
.checkDendroClusterFormat<-function(dendro,checkLabels=TRUE){
    data.cl<-phylobase::tdata(dendro)
    if(!all(.clusterDendroColumns %in% names(data.cl) )){
        return("dendro_clusters must have data with column names:",paste(.clusterDendroColumns,sep=","))
    }
    tp<-phylobase::tipLabels(dendro)
    if(checkLabels && any(sort(as.numeric(gsub("T","",tp)))!=sort(seq_along(tp))))
        return("dendro_clusters cannot have labels for the tips; user-defined labels for the tips (i.e. clusters) should be stored in the clusterLegend")
    if(any(is.na(data.cl$Position))) 
        return("dendro_clusters cannot have NA values in Position variable")
    if( any(is.na(data.cl$NodeId))) 
        return("dendro_clusters cannot have NA values in Node Id variable")
    return(TRUE)
}
.checkDendroSamplesFormat<-function(dendro,checkLabels=TRUE){
    data.cl<-phylobase::tdata(dendro,type="all")
    all(names(data.cl)%in% .clusterSampleColumns)
    if(!all(.clusterSampleColumns %in% names(data.cl) )){
        return("dendro_samples must have data with column names:",paste(.clusterDendroColumns,sep=","))
    }
    if(any(is.na(data.cl$Position))) 
        return("dendro_samples cannot have NA values in Position variable")
    if(checkLabels && any(!is.na(phylobase::nodeLabels(dendro)))) 
        return("dendro_samples cannot have labels for the nodes; the labels for the nodes should be in the cluster dendrogram")
    if(checkLabels){
        tp<-phylobase::tipLabels(dendro)
        if(any(sort(as.numeric(gsub("T","",tp)))!=sort(seq_along(tp)))) return("dendro_samples cannot have labels for the tips; the labels for the tips should be the colnames of the object")
    } 
    data.cl<-phylobase::tdata(dendro,type="tip")
    if(any(is.na(data.cl$SampleIndex))) 
        return("dendro_samples cannot have NA values for tips in SampleIndex variable")
    return(TRUE)
    
}
.checkDendrogram<-function(object){
    ############
    ##Check dendrogram
    ############
    if(!is.null(object@dendro_clusters)){
        if(is.na(dendroClusterIndex(object))) return("if dendrogram slots are filled, must have corresponding dendro_index defined.")
        dcluster<-clusterMatrix(object)[,dendroClusterIndex(object)]
        if(phylobase::nTips(object@dendro_clusters) != max(dcluster)) {
            return("dendro_clusters must have the same number of leaves as the number of (non-negative) clusters")
        }
        ch<-.checkDendroClusterFormat(object@dendro_clusters)
        if(!is.logical(ch)) return(ch)
        
        #further checks that require full CE object:
        data.cl<-phylobase::tdata(object@dendro_clusters,type="all")
        if(any(!gsub("ClusterId","",na.omit(data.cl$ClusterIdDendro)) %in% as.character(object@clusterMatrix[,object@dendro_index]))) return("ClusterIdDendro information in dendrogram slot must match the corresponding cluster ids in clustering defined by dendro_index slot")
        if(any(!gsub("ClusterId","",na.omit(data.cl$ClusterIdMerge)) %in% as.character(object@clusterMatrix[,object@merge_index]))) return("ClusterIdMerge information in dendrogram slot must match the corresponding cluster ids in clustering defined by merge_index slot")
    }
    else{
        if(!is.null(object@dendro_samples)) return("dendro_clusters should not be null if dendro_samples is non-null") #if comment out now optional to have samples
    }
    if(!is.null(object@dendro_samples)){
        if(phylobase::nTips(object@dendro_samples) != NCOL(object)) {
            return("dendro_samples must have the same number of leaves as the number of samples")
        }
        ch<-.checkDendroSamplesFormat(object@dendro_samples)
        if(!is.logical(ch)) return(ch)
        
    }
    else{
        if(!is.null(object@dendro_clusters)) return("dendro_samples should not be null if dendro_clusters is non-null") #if commented out, makes it optional to have samples
    }
    
    return(TRUE)
}


.checkMerge<-function(object){
  ############
  ##Check merge related slots
  ############
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
    if(!baseMergeMethod %in% .validMergeMethods) return(paste("merge_method must be one of available merge methods:", paste(.validMergeMethods,collapse=",")," (with possibility of fold-change added to method name for 'adjP')"))
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
    if(!all(requireColumns %in% namesToCheck))
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
