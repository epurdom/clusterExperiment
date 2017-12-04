#####Function to calculate the reducedDim matrices.


#' @rdname reduceFunctions
#' @param reducedDims a vector of character values indicating the methods of dimensionality reduction to be performed. Currently only "PCA" is implemented.
#' @param maxDims Numeric vector of integer giving the number of PC dimensions to calculate. 
#'   \code{maxDims} can also take values between (0,1) to indicate keeping the
#'   number of dimensions necessary to account for that proportion of the
#'   variance. \code{maxDims} should be of same length as \code{reducedDims}, indicating 
#'   the number of dimensions to keep for each method (if \code{maxDims} is of length 1, 
#'   the same number of dimensions will be used for each). 
#' @param ... Values passed on the the 'SingleCellExperiment' method.
#' @inheritParams transformData
#' @details The PCA method uses either \code{prcomp} from the \code{stats} package or  \code{svds} from the \code{RSpectra} package to perform PCA. Both are called on \code{t(assay(x))} with \code{center=TRUE} and \code{scale=TRUE} (i.e. the feature are centered and scaled), so that
#'   it is performing PCA on the correlation matrix of the features. 
#'
#' @return \code{makeReducedDims} returns a \code{\link{SingleCellExperiment}} containing the calculated dimensionality reduction in the \code{reduceDims} with names corresponding to the name given in \code{reducedDims}.
#' @examples
#' data(simData)
#' listBuiltInReducedDims()
#' scf<-makeReducedDims(simData, reducedDims="PCA", maxDims=3)
#' scf
#' @export
#' @aliases makeReducedDims,SingleCellExperiment-method makeReducedDims
#' @importFrom matrixStats rowVars
setMethod(
  f = "makeReducedDims",
  signature = "SingleCellExperiment",
  definition = function(object,reducedDims="PCA",maxDims=500,transFun=NULL,isCount=FALSE)
{

  ###################
  ##Check user inputs
  ###################
  #check valid options for reducedDims
  validDim<-listBuiltInReducedDims()
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
#' @rdname reduceFunctions
#' @export
setMethod(
  f = "makeReducedDims",
  signature = "matrix",
  definition = function(object,...)
{
	makeReducedDims(SummarizedExperiment(object),...)
}
)
#' @rdname reduceFunctions
#' @export
setMethod(
  f = "makeReducedDims",
  signature = "SummarizedExperiment",
  definition = function(object,...)
{
	makeReducedDims(as(object,"SingleCellExperiment"),...)
}
)
#' @rdname reduceFunctions
#' @export
setMethod(
  f = "makeReducedDims",
  signature = "ClusterExperiment",
  definition = function(object,...)
{
	if(any(c("transFun","isCount") %in% names(list(...)))) 
		stop("The internally saved transformation function of a ClusterExperiment object must be used when given as input and setting 'transFun' or 'isCount' for a 'ClusterExperiment' is not allowed.")  
	out<-makeReducedDims(as(object,"SingleCellExperiment"),transFun=transformation(object),...)
	return(.addBackSEInfo(newObj=object,oldObj=out))
}
)

#' @rdname reduceFunctions
#' @export
listBuiltInReducedDims<-function(){c("PCA")}

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



