#' @name transformData
#' @title Transform the original data in a ClusterExperiment object
#'
#' @description Provides the transformed data
#'
#' @param object a matrix, SummarizedExperiment, SingleCellExperiment or ClusterExperiment object.
#' @param transFun a transformation function to be applied to the data. If the transformation 
#' applied to the data creates an error or NA values, then the function will throw an error.
#' If object is of class \code{ClusterExperiment}, the stored transformation will be used 
#' and giving this parameter will result in an error.
#' @param isCount if \code{transFun=NULL}, then \code{isCount=TRUE} will determine the
# transformation as defined by \code{function(x){log2(x+1)}}, and \code{isCount=FALSE} 
#' will give a transformation function \code{function(x){x}}. Ignored if \code{transFun=NULL}.
#' If object is of class \code{ClusterExperiment}, the stored transformation will be used 
#' and giving this parameter will result in an error.
#' @param ... Values passed on the the 'matrix' method.
#' @details The data matrix defined by \code{assay(x)} is transformed based on
#'   the transformation function either defined in x (in the case of a 
#' \code{ClusterExperiment} object) or by user given values for other classes. 
#'
#'
#' @examples
#' mat <- matrix(data=rnorm(200), ncol=10)
#' mat[1,1] <- -1 #force a negative value
#' labels <- gl(5, 2)
#' cc <- ClusterExperiment(mat, as.numeric(labels), transformation =
#' function(x){x^2}) #define transformation as x^2
#' z<-transformData(cc) 
#' @aliases transformData,matrix-method
#' @export
setMethod(
  f = "transformData",
  signature = "matrix",
  definition = function(object,transFun=NULL,isCount=FALSE) {
	  transFun<-.makeTransFun(transFun=transFun,isCount=isCount)
	  x <- try(transFun(object), silent=TRUE)
	  if(inherits(x, "try-error"))
	    stop("User-supplied `transFun` produces the following error on the input data matrix:\n",x)
	  if(anyNA(x))
	    stop("User-supplied `transFun` produces NA values")
	  return(x)
  }
)
#' @export
#' @rdname transformData
setMethod(
  f = "transformData",
  signature = "ClusterExperiment",
  definition = function(object,...) {
  	if(any(c("transFun","isCount") %in% names(list(...)))) 
  		stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")  
	  return(transformData(assay(object),transFun=transformation(object)))
  }
)
#' @export
#' @rdname transformData
setMethod(
  f = "transformData",
  signature = "SingleCellExperiment",
  definition = function(object,...) {
	  return(transformData(assay(object),...))
  }
)
#' @export
#' @rdname transformData
setMethod(
  f = "transformData",
  signature = "SummarizedExperiment",
  definition = function(object,...) {
	  return(transformData(as(object,"SingleCellExperiment"),...))
  }
)

#small function to uniformally return transformation function from combination of transFun and isCount
.makeTransFun<-function(transFun=NULL,isCount=FALSE){
  if(is.null(transFun)){
    transFun <- if(isCount) function(x){log2(x+1)} else function(x){x}
  }
  return(transFun)
}


#####Function to calculate the reducedDim matrices.


#' @name makeDimReduce
#' @title Calculate dimensionality reduction of data
#' @description A Function for performing and storing common dimensionality reduction techniques
#' @param reducedDims a vector of character values indicating the methods of dimensionality reduction to be performed. Currently only "PCA" is implemented.
#' @param object input to use for the data for dimensionality reduction. Can be matrix, SummarizedExperiment, SingleCellExperiment, or ClusterExperiment object
#' @param maxDims Numeric vector of integer giving the number of PC dimensions to calculate. 
#'   \code{maxDims} can also take values between (0,1) to indicate keeping the
#'   number of PCA dimensions necessary to account for that proportion of the
#'   variance. \code{maxDims} should be of same length as \code{reducedDims}, indicating the number of dimensions to keep for each method (if \code{maxDims} is of length 1, the same number of dimensions will be used for each). 
#' @param ... Values passed on the the 'SingleCellExperiment' method.
#' @inheritParams transformData
#' @return a SingleCellExperiment object with the indicated diminsionality reduction methods stored in the \code{reduceDims} slot.
#' @details The PCA method uses either \code{prcomp} from the \code{stats} package or  \code{svds} from the \code{RSpectra} package to perform PCA. Both are called on \code{t(assay(x))} with \code{center=TRUE} and \code{scale=TRUE} (i.e. the feature are centered and scaled), so that
#'   it is performing PCA on the correlation matrix of the features. 
#'
#' @return A \code{\link{SingleCellExperiment}} containing the dimensionality reduction in the corresponding slots with names corresponding to the name given in \code{reducedDims}.
#' @examples
#' data(simData)
#' listBuiltInDimReduce()
#' scf<-makeDimReduce(simData, reducedDims="PCA", maxDims=3)
#' scf
#' @export
#' @aliases makeDimReduce,SingleCellExperiment-method
#' @importFrom matrixStats rowVars
setMethod(
  f = "makeDimReduce",
  signature = "SingleCellExperiment",
  definition = function(object,reducedDims="PCA",maxDims=500,transFun=NULL,isCount=FALSE)
{

  ###################
  ##Check user inputs
  ###################
  #check valid options for reducedDims
  validDim<-listBuiltInDimReduce()
  reducedDims<-unique(reducedDims)
  if(length(maxDims)==1) maxDims<-rep(maxDims,length=length(reducedDims))
  if(length(maxDims)!=length(reducedDims)) stop("'maxDims' must be of same length as 'reducedDims'")
	  
  ######Check dimensions and valid argument
  for(dr in reducedDims){
	  dr<-match.arg(dr,validDim) 
	  if(is.na(maxDims) || maxDims>NROW(object) || maxDims > NCOL(object)){
		  if(!is.na(maxDims) & (maxDims>NROW(object) || maxDims > NCOL(object)))
			  warning("User requested more dimensionality reduction dimensions than the minimimum of number of rows and columns. Will return all dimensions.")
		  maxDims<-min(c(NROW(object),NCOL(object)))
	  }
	  if(maxDims<=0)  stop("the number of reducedDims dimensions must be a value strictly greater than 0")


  }
  ###################
  ##Clean up data:
  ###################
  #transform data
  x<-transformData(object,transFun=transFun,isCount=isCount)
  #---------
  #Check zero variance genes before doing reducedDims:
  #---------
  rowvars <- matrixStats::rowVars(x)
  if(any(rowvars==0)) {
    if(all(rowvars==0)) {
      stop("All features have zero variance.")
    }
    warning("Found features with zero variance.\nMost likely these are features with 0 across all samples.\nThey will be removed from dimensionality reduction step.")
  }

  ###################
  ##Do loop over reducedDims values:
  ###################
  currErrors<-c()
  for(kk in 1:length(reducedDims)){
	  dr<-reducedDims[[kk]]
	  md<-maxDims[[kk]]
	  isPct <- md < 1
	  #check if already calculated
	  #note, currently no way to check if have already done if md<1
	  if(dr %in% reducedDimNames(object)){
		  if(!isPct & md<=ncol(reducedDim(object,type=dr))) next
	  }
	  #-------------
	  # if add other functions, add if statements here
	  if(dr=="PCA") out<-try(.pcaDimRed(x,md=md,isPct=isPct,rowvars=rowvars))
	  ##-------
	  	  
	  if(!inherits(out,"try-error")) reducedDim(object,reducedDims) <- out
	  else{
		  currErrors<-c(currErrors,paste("\t",dr,":",out,sep=""))
	  }	  
  }
  if(length(currErrors)>0){
	  if(length(currErrors)==length(reducedDims)) 
		  stop(paste("No dimensionality reduction techniques were successful:",currErrors,sep="\n"))
	  else{
	  	warning(paste("The following dimensionality reduction techniques hit errors:",currErrors,sep="\n"))
	  }
  }
  return(object)

}
)
#' @rdname makeDimReduce
#' @export
setMethod(
  f = "makeDimReduce",
  signature = "matrix",
  definition = function(object,...)
{
	makeDimReduce(SummarizedExperiment(object),...)
}
)
#' @rdname makeDimReduce
#' @export
setMethod(
  f = "makeDimReduce",
  signature = "SummarizedExperiment",
  definition = function(object,...)
{
	makeDimReduce(as(object,"SingleCellExperiment"),...)
}
)
#' @rdname makeDimReduce
#' @export
setMethod(
  f = "makeDimReduce",
  signature = "ClusterExperiment",
  definition = function(object,...)
{
	if(any(c("transFun","isCount") %in% names(list(...)))) 
		stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")  
	out<-makeDimReduce(as(object,"SingleCellExperiment"),transFun=transformation(object),...)
	return(.addBackSEInfo(newObj=object,oldObj=out))
}
)

#' @rdname makeDimReduce
#' @export
listBuiltInDimReduce<-function(){c("PCA")}

#' @importFrom RSpectra svds
.pca <- function(x, center=TRUE, scale=FALSE, k) {
  svd_raw <- svds(scale(x, center=center, scale=scale), k=k, nu=k, nv=0)
  pc_raw <- svd_raw$u %*% diag(svd_raw$d, nrow = length(svd_raw$d))
  rownames(pc_raw) <- rownames(x)
  return(pc_raw)
}
#' @importFrom stats prcomp
.pcaDimRed<-function(x,md,isPct,rowvars){	
	if(isPct) {
		prcObj<-stats::prcomp(t(x[which(rowvars>0),]),center=TRUE,scale=TRUE)
		prvar<-prcObj$sdev^2 #variance of each component
		prvar<-prvar/sum(prvar)
		prc<-prcObj$x
		if(NROW(prc) != NCOL(x)) stop("Internal error in coding of principal components.")
		md <- which(cumsum(prvar)>md)[1] #pick first pca coordinate with variance > value
		prc <- prc[,seq_len(md)]
	}
	else {
		prc <- .pca(t(x[which(rowvars>0),]), center=TRUE, scale=TRUE, k=md)
		if(any(md > NROW(prc)))
		  stop("Internal error in coding of principal components.")
	}
	colnames(prc)<-paste("PC",1:ncol(prc),sep="") #make them match prcomp
	return(prc)
}







	
#' @name makeFilterStats
#' @title Calculate filtering statistics
#' @description Function for calculating, per row (gene), built-in statistical 
#'    functions that might be used for filtering. 
#' @param object object from which user wants to calculate per-row statistics
#' @param filterStats character vector of statistics to calculate. 
#' 	  Must be one of the character values given by \code{listBuildInFilterStats()}.
#' @param ... Values passed on the the 'SingleCellFilter' method.
#' @inheritParams makeDimReduce
#' @inheritParams transformData
#'
#' @examples
#' data(simData)
#' listBuiltInFilterStats()
#' scf<-makeFilterStats(simData,filterStats=c("var","mad"))
#' scf
#' scfFiltered<-filterData(scf,type="mad",percentile=10)
#' scfFiltered
#' assay(scfFiltered)[1:10,1:10]
#' @rdname makeFilterStats
#' @aliases makeFilterStats,SingleCellFilter-method
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "SingleCellFilter",
  definition = function(object,filterStats=listBuiltInFilterStats(),transFun=NULL,isCount=FALSE)
{

  ###################
  ##Check user inputs
  ###################
  #check valid options for reducedDims
  validDim<-listBuiltInFilterStats()
  filterStats<-unique(filterStats)
  if(!all(filterStats %in% validDim)){
	  stop("Not all of the filterStats given are valid")
  }
  
  ###################
  ##Clean up data:
  ###################
  #transform data
  x<-transformData(object,transFun=transFun,isCount=isCount)

  ###################
  ##Do loop over filterStats values:
  ###################
  if('abscv' %in% filterStats){
	  #origfilterStats<-filterStats
  #if abscv, add var and mean since calculating anyway and make abscv last one
	  filterStats<-unique(c(filterStats,"var","mean"))
	  doCV<-TRUE
	  whCV<-grep('abscv',filterStats)
	  filterStats<-c(filterStats[-whCV])
  }
  else doCV<-FALSE
  filterStatData<-sapply(filterStats,function(statName){
	  f<-.matchToStats[[statName]]
	  f(x)
  })
  if(doCV){
	  filterStatData<-cbind(filterStatData, "abscv"=sqrt(filterStatData[,"var"])/abs(filterStatData[,"mean"]))
	  #filterStatData<-filterStatData[,origfilterStats] #put it in order, though user shouldn't depend on it.
  }
  filterStats(object)<-filterStatData #should leave in place existing ones, update conflicting ones, and add new ones!
  return(object)

}
)
#' @rdname makeFilterStats
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "matrix",
  definition = function(object,...)
{
	makeFilterStats(SummarizedExperiment(object),...)
}
)
#' @rdname makeFilterStats
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "SummarizedExperiment",
  definition = function(object,...)
{
	makeFilterStats(as(object,"SingleCellExperiment"),...)
}
)
#' @rdname makeFilterStats
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "SingleCellExperiment",
  definition = function(object,...)
{
	makeFilterStats(as(object,"SingleCellFilter"),...)
}
)
#' @rdname makeFilterStats
#' @export
#' @param whichClusterIgnoreUnassigned indicates clustering that should be used to filter out unassigned samples from the calculations. If \code{NULL} no filtering of samples will be done. See details for more information.
#' @details \code{whichClusterIgnoreUnassigned} is only an option when applied to a
#' \code{ClusterExperiment} classs and indicates that the filtering statistics should be
#' calculated based on samples that are unassigned by the
#'  designated clustering. The name given to the filter in this case is of the form
#'  \code{<filterStats>_<clusterLabel>}, i.e. the clustering label of the clustering is
#'  appended to the standard name for the filtering statistic.
#'

setMethod(
  f = "makeFilterStats",
  signature = "ClusterExperiment",
  definition = function(object,whichClusterIgnoreUnassigned=NULL,filterStats=listBuiltInFilterStats(),...)
{
	if(any(c("transFun","isCount") %in% names(list(...)))) 
		stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")  
	filterStats<-unique(filterStats)
	if(!is.null(whichClusterIgnoreUnassigned)){
		whCluster<-.TypeIntoIndices(object,whichClusterIgnoreUnassigned)
		if(length(whCluster)>1) warning("'whichClusterIgnoreUnassigned' corresponds to multiple clusterings. Ignoring input")
		else if(length(whCluster)==0) warning("'whichClusterIgnoreUnassigned' does not correspond to a clustering. Ignoring input")
		else{
			#give new names to filters to indicate based on clustering.
			newNames<-paste(filterStats,clusterLabels(object)[whCluster],sep="_")
			whDo<-which(!newNames %in% filterNames(object))
			if(length(whDo)>0){
				whAssigned<-which(clusterMatrix(object)[,whCluster]>0)
				if(length(whAssigned)>0){
					out<-makeFilterStats(object[,whAssigned],filterStats=filterStats[whDo],...)
					whNew<-match(filterStats[whDo],filterNames(out))
					filterNames(out)[whNew]<-newNames[whDo]	
				}
				else 
					stop("All samples are unassigned for clustering", clusterLabels(object)[whCluster])
				
			}
		}
	}
	else{
		out<-makeFilterStats(as(object,"SingleCellFilter"),filterStats=filterStats,transFun=transformation(object),...)
	}
#	
	filterStats(object)<-filterStats(out)
	return(object)
}
)
#' @rdname makeFilterStats
#' @export
listBuiltInFilterStats<-function(){c('var', 'abscv', 'mad','mean','iqr','median')}
#' @importFrom matrixStats rowVars rowMeans2 rowMads rowMedians rowIQRs
.matchToStats<-SimpleList(
	'var'=matrixStats::rowVars,
	'mad'=matrixStats::rowMads,
	'mean'=matrixStats::rowMeans2,
	'iqr'=matrixStats::rowIQRs,
	'median'=matrixStats::rowMedians)


#' @rdname makeFilterStats
#' @aliases filterData
#' @param type The type of filtering statistic to use to filter. 
#' @param cutoff numeric. A value at which to filter the rows (genes) for the test statistic
#' @param percentile numeric. Either a number between 0,1 indicating what percentage of the rows (genes) to keep or an integer value indicated the number of rows (genes) to keep
#' @param absolute whether to take the absolute value of the filter statistic
#' @param keepLarge logical whether to keep rows (genes) with large values of the test statistic or small values of the test statistic. 
#' @details Note that \code{filterData} returns a SingleCellFilter object. To get the actual data out use either assay or \code{\link{transformData}} if transformed data is desired.
#' @return A SingleCellFilter object with the rows (genes) removed based on filters
#' @export
#' @importFrom stats quantile
setMethod( "filterData","SingleCellFilter",
	function(object,type,cutoff,percentile, absolute=FALSE,keepLarge=TRUE){
	stat<-if(absolute) abs(filterStats(object,type)) else filterStats(object,type)
	if(missing(cutoff) & missing(percentile)) stop("must provide one of cutoff or percentile")
	if(!missing(cutoff) & !missing(percentile)) stop("can only provide one of cutoff or percentile")
	if(!missing(cutoff)){
		whKeep<- if(keepLarge) which(stat>cutoff) else which(cutoff > stat)
	}
	if(!missing(percentile)){
		if(0<percentile & percentile <1){
			quantile<- quantile(stat,probs=if(keepLarge) percentile else 1-percentile)
			whKeep<-if(keepLarge) which(stat>quantile) else which(stat<quantile)
		}
		else{
			if(is.na(percentile) || percentile>=1){
				if(is.na(percentile) || percentile>NROW(object)){
					warning("the number of most features requested after filtering is either missing or larger than the number of features. Will not do any filtering")
					whKeep<-1:NROW(object)
				}
				else whKeep<- order(stat,decreasing=ifelse(keepLarge,TRUE,FALSE))[1:percentile]
			}
			else stop("Invalid value for percentile. Must be either between 0,1 or a positive integer number to keep")
		}
	}
	object[whKeep,]

})	