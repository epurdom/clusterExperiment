############
## Checks that make up validity of clusterExperiment:
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
  if(any(is.na(assay(object)))) {
    return("NA values are not allowed.")
  }
  return(TRUE)
}
.checkTransform<-function(object){
    #############
    #check transform
    #############
    tX <- try(transform(object),silent=TRUE)
    if(inherits(tX, "try-error")){
      stop(paste("User-supplied `transformation` produces error on the input data
                 matrix:\n",tX))
    }
    if(any(is.na(tX))) {
      return("NA values after transforming data matrix are not allowed.")
    }
	return(TRUE)
}

.checkClusterMatrix<-function(object){
    ############
    #Check clusterMatrix
    ############

    if(!all(is.na((object@clusterMatrix))) &
       !(NROW(object@clusterMatrix) == NCOL(object))) {
      return("If present, `clusterMatrix` must have as many row as cells.")
    }
    if(!is.numeric(object@clusterMatrix)) {
      return("`clusterMatrix` must be a numeric matrix.")
    }

    if(NCOL(object@clusterMatrix)!= length(object@clusterTypes)) {
      return("length of clusterTypes must be same as NCOL of the clusterMatrix")
    }

    if(NCOL(object@clusterMatrix)!= length(object@clusterInfo)) {
      return("length of clusterInfo must be same as NCOL of the clusterMatrix")
    }
        #check internally stored as integers
        testConsecIntegers<-apply(object@clusterMatrix,2,function(x){
          whCl<-which(!x %in% c(-1,-2))
          uniqVals<-unique(x[whCl])
          return(all(sort(uniqVals)==1:length(uniqVals)))
        })
     if(!all(testConsecIntegers)) return("the cluster ids in clusterMatrix must be stored internally as consecutive integer values")
	 
	 return(TRUE)
}
.checkClusterLabels<-function(object){
	 if(!all(is.na(object@clusterMatrix))){
	     if(is.null(colnames(object@clusterMatrix))) return("clusterMatrix must have column names")
	     if(any(duplicated(colnames(object@clusterMatrix)))) return("clusterMatrix must have unique column names")
 	
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
      if(!is.null(object@dendro_clusters)) return("dendro_samples should not be null if dendro_clusters is non-null")
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
      if(!is.null(object@dendro_samples)) return("dendro_clusters should not be null if dendro_samples is non-null")
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

.checkClusterLegend<-function(object){
    ####
    #test that @clusterLegend is proper form
    ####
	if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
	    if(!is.null(names(object@clusterLegend))) return("clusterLegend should not have names")
	    if(length(object@clusterLegend) != NCOL(object@clusterMatrix)) {
	      return("`clusterLegend` must be list of same length as NCOL of
	               `clusterMatrix`")
	    }
	    testIsMatrix <- sapply(object@clusterLegend,
	                           function(x) {!is.null(dim(x))})
	    if(!all(testIsMatrix)) {
	      return("Each element of `clusterLegend` list must be a matrix")
	    }
	    testColorRows <- sapply(object@clusterLegend, function(x){nrow(x)})
	    testClusterMat <- apply(object@clusterMatrix, 2, function(x) {length(unique(x))})
	    if(!all(testColorRows == testClusterMat)) {
	      return("each element of `clusterLegend` must be matrix with number of
	               rows equal to the number of clusters (including -1 or -2 values)
	               in `clusterMatrix`")
	    }
	    testColorCols1 <- sapply(object@clusterLegend, function(x) {
	      "color" %in% colnames(x)})
	    testColorCols2 <- sapply(object@clusterLegend, function(x) {
	      "clusterIds" %in% colnames(x)})
	    testColorCols3 <- sapply(object@clusterLegend, function(x) {
	      "name" %in% colnames(x)})
	    if(!all(testColorCols1) || !all(testColorCols2) || !all(testColorCols3)) {
	      return("each element of `clusterLegend` must be matrix with at least 3
	             columns, and at least 3 columns have names `clusterIds`,
	             `color` and `name`")
	    }
	#     testUniqueName <- sapply(object@clusterLegend, function(x) {
	#       any(duplicated(x[,"name"]))})
	#     if(any(testUniqueName)) return("the column")
	    testColorCols1 <- sapply(object@clusterLegend, function(x){is.character(x)})
	    if(!all(testColorCols1)) {
	      return("each element of `clusterLegend` must be matrix of character
	             values")
	    }
	    testColorCols1 <- sapply(1:length(object@clusterLegend), function(ii){
	      col<-object@clusterLegend[[ii]]
	      x<-object@clusterMatrix[,ii]
	      y<-as.numeric(col[,"clusterIds"])
	      all(y %in% x)
	    })
	    if(!all(testColorCols1)) {
	      return("each element of `clusterLegend` must be matrix with column
	             `clusterIds` matching the corresponding integer valued
	             clusterMatrix values")
	    }
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
	if(any(!object@orderSamples %in% 1:NCOL(assay(object)))) {
	  return("`orderSamples` must be values between 1 and the number of samples.")
	}
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
#' myCF<-clusterFunction(clusterFUN=goodFUN, inputType="X",algorithmType="K", outputType="vector")
#' badFUN<-function(x,diss,k,checkArgs,cluster.only,...){cluster::pam(x=x,k=k)}
#' internalFunctionCheck(badFUN,inputType="X",algorithmType="K",outputType="vector")
#' @details \code{internalFunctionCheck} is the function that is called by the 
#'   validity check of the \code{clusterFunction} constructor (if 
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