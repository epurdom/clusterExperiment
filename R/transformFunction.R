#' @name transformData
#' @title Transform the original data in a ClusterExperiment object
#'
#' Provides the transformed data (as defined by the object)
#'
#' @param _data a matrix, SummarizedExperiment, SingleCellExperiment or ClusterExperiment object.
#'
#' @details The data matrix defined by \code{assay(x)} is transformed based on
#'   the transformation function defined in x. 
#'
#'
#' @examples
#' mat <- matrix(data=rnorm(200), ncol=10)
#' mat[1,1] <- -1 #force a negative value
#' labels <- gl(5, 2)
#' cc <- ClusterExperiment(mat, as.numeric(labels), transformation =
#' function(x){x^2}) #define transformation as x^2
#' z<-transformData(cc) 
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
  definition = function(object) {
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


#####Function to calculate the dimReduce matrices.


#' @name makeDimReduce
#' @title Calculate dimensionality reduction of data
#' @description A Function for performing and storing common dimensionality reduction techniques
#' @param dimReduce a vector of character values indicating the methods of dimensionality reduction to be performed. Currently only "PCA" is implemented.
#' @param object input to use for the data for dimensionality reduction. Can be matrix, SummarizedExperiment, SingleCellExperiment, or ClusterExperiment object
#' @param maxDims Numeric vector of integer giving the number of PC dimensions to calculate. 
#'   \code{maxDims} can also take values between (0,1) to indicate keeping the
#'   number of PCA dimensions necessary to account for that proportion of the
#'   variance. \code{maxDims} should be of same length as \code{dimReduce}, indicating the number of dimensions to keep for each method (if \code{maxDims} is of length 1, the same number of dimensions will be used for each). 
#' @param transFun a function that will be used to transform the data before performing dimensionality reduction. If \code{object} is a \code{ClusterExperiment} object, then the value of \code{transformation} slot will be used and the user cannot give a value for \code{transFun}.
#' @param isCount logical. If \code{transFun} is not given, and \code{isCount=TRUE}, then the \code{transFun} will be assumed to be \code{log(x+1)}.
#' @return a SingleCellExperiment object with the indicated diminsionality reduction methods stored in the \code{reduceDims} slot.
#' @details The PCA method uses either \code{prcomp} from the \code{stats} package or  \code{svds} from the \code{RSpectra} package to perform PCA. Both are called on \code{t(assay(x))} with \code{center=TRUE} and \code{scale=TRUE} (i.e. the feature are centered and scaled), so that
#'   it is performing PCA on the correlation matrix of the features. 
#'
#' @return A \code{\link{SingleCellExperiment}} containing the dimensionality reduction in the corresponding slots with names corresponding to the name given in \code{dimReduce}.
#' @examples
#' mat <- matrix(data=rnorm(200), ncol=10)
#' mat[1,1] <- -1 #force a negative value
#' labels <- gl(5, 2)
#' cc <- ClusterExperiment(mat, as.numeric(labels), transformation =
#' function(x){x^2}) #define transformation as x^2
#' #will transform data based on saved transformation
#' x <- makeDimReduce(cc, dimReduce="PCA", nPCADims=3)
#' @export
#' @importFrom matrixStats rowVars
setMethod(
  f = "makeDimReduce",
  signature = "SingleCellExperiment",
  definition = function(object,dimReduce="PCA",maxDims=500,transFun=NULL,isCount=FALSE)
{

  ###################
  ##Check user inputs
  ###################
  #check valid options for dimReduce
  validDim<-listBuiltInDimReduce()
  dimReduce<-unique(dimReduce)
  if(length(maxDims)==1) maxDims<-rep(maxDims,length=length(dimReduce))
  if(length(maxDims)!=length(dimReduce)) stop("'maxDims' must be of same length as 'dimReduce'")
	  
  ######Check dimensions and valid argument
  for(dr in dimReduce){
	  dr<-match.arg(dr,validDim) 
	  if(is.na(maxDims) || maxDims>NROW(object) || maxDims > NCOL(object)){
		  maxDims<-min(c(NROW(object),NCOL(object)))
	  }
	  if(maxDims<=0)  stop("the number of dimReduce dimensions must be a value strictly greater than 0")


  }
  ###################
  ##Clean up data:
  ###################
  #transform data
  x<-transformData(object,transFun=transFun,isCount=isCount)
  #---------
  #Check zero variance genes before doing dimReduce:
  #---------
  rowvars <- matrixStats::rowVars(x)
  if(any(rowvars==0)) {
    if(all(rowvars==0)) {
      stop("All features have zero variance.")
    }
    warning("Found features with zero variance.\nMost likely these are features with 0 across all samples.\nThey will be removed from dimensionality reduction step.")
  }

  ###################
  ##Do loop over dimReduce values:
  ###################
  currErrors<-c()
  for(kk in 1:length(dimReduce)){
	  dr<-dimReduce[[kk]]
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
	  	  
	  if(!inherits(out,"try-error")) reducedDim(object,dimReduce) <- out
	  else{
		  currErrors<-c(currErrors,paste("\t",dr,":",out,sep=""))
	  }	  
  }
  if(length(currErrors)>0){
	  if(length(currErrors)==length(dimReduce)) 
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
	out<-makeDimReduce(as(object,"SingleCellExperiment"),transFun=transformation(object),...)
	return(.addBackSEInfo(newObj=object,oldObj=out))
}
)

#' @rdname makeDimReduce
#' @export
listBuiltInDimReduce<-function(){c("PCA")}

#' @importFrom RSpectra svds
#' @importFrom stats prcomp
.pcaDimRed<-function(x,md,isPct,rowvars){
	.pca <- function(x, center=TRUE, scale=FALSE, k) {
	  svd_raw <- svds(scale(x, center=center, scale=scale), k=k, nu=k, nv=0)
	  pc_raw <- svd_raw$u %*% diag(svd_raw$d, nrow = length(svd_raw$d))
	  rownames(pc_raw) <- rownames(x)
	  return(pc_raw)
	}
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







#' @rdname makeFilterStats
#' @export
listBuiltInFilterStats<-function(){c('var', 'cv', 'mad','mean','iqr','median')}
#' @importFrom matrixStats rowVars rowMeans2 rowMads rowMedians rowIQRs
.matchToStats<-SimpleList(
	'var'=matrixStats::rowVars,
	'mad'=matrixStats::rowMads,
	'mean'=matrixStats::rowMeans2,
	'iqr'=matrixStats::rowIQRs,
	'median'=matrixStats::rowMedians)
	
#' @export
setMethod(
  f = "makeFilterStats",
  signature = "SingleCellExperiment",
  definition = function(object,filterStats=listBuiltInFilterStats(),transFun=NULL,isCount=FALSE)
{

  ###################
  ##Check user inputs
  ###################
  #check valid options for dimReduce
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
  currErrors<-c()
  if('cv' %in% filterStats){
	  #if cv, add var and mean since calculating anyway and make cv last one
	  filterStats<-unique(c(filterStats,"var","mean"))
	  doCV<-TRUE
	  whCV<-grep('cv',filterStats)
	  filterStats<-c(filterStats[-whCV])
  }
  else doCV<-FALSE
  filterStatData<-sapply(filterStats,function(statName){
	  f<-.matchToStats[[statName]]
	  f(x)
  })
  if(doCV){
	  filterStatData<-cbind(filterStatData, "cv"=sqrt(filterStatData[,"var"])/filterStatData[,"mean"])
  }
  
  # for(st in filterStats){
  # 	  if(st!="cv"){
  # 		  function
  # 	  }
  # 	  #-------------
  # 	  # if add other functions, add if statements here
  # 	  if(dr=="PCA") out<-try(.pcaDimRed(x,md=md,isPct=isPct,rowvars=rowvars))
  # 	  ##-------
  #
  # 	  if(!inherits(out,"try-error")) reducedDim(object,dimReduce) <- out
  # 	  else{
  # 		  currErrors<-c(currErrors,paste("\t",dr,":",out,sep=""))
  # 	  }
  # }
  # if(length(currErrors)>0){
  # 	  if(length(currErrors)==length(dimReduce))
  # 		  stop(paste("No dimensionality reduction techniques were successful:",currErrors,sep="\n"))
  # 	  else{
  # 	  	warning(paste("The following dimensionality reduction techniques hit errors:",currErrors,sep="\n"))
  # 	  }
  # }
  return(object)

}
)
#' @rdname makeDimReduce
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
  signature = "ClusterExperiment",
  definition = function(object,...)
{
	out<-makeFilterStats(as(object,"SingleCellExperiment"),transFun=transformation(object),...)
	return(.addBackSEInfo(newObj=object,oldObj=out))
}
)


#' @details \code{ignoreUnassignedVar} has no impact for PCA reduction, which
#'   will always use all samples. At all times, regardless of the value of
#'   \code{ignoreUnassignedVar}, a matrix with the same number of columns of
#'   \code{assay(x)} (i.e. the same number of samples) will be returned.
#' @details  \code{dimReduce}, \code{nPCADims}, \code{nVarDims} can all be a
#'   vector of values, in which case a list will be returned with the
#'   appropriate datasets as elements of the list.
#' @param nVarDims Numeric (integer) vector giving the number of features (e.g.
#'   genes) to keep, based on variance/cv/mad variability.
#' @param dimReduce Character vector specifying the dimensionality reduction to
#'   perform, any combination of 'none', 'PCA', 'var', 'cv', and 'mad'. See details.
#' @param ignoreUnassignedVar logical indicating whether dimensionality reduction
#'   via top feature variability (i.e. 'var','cv','mad') should ignore
#'   unassigned samples in the primary clustering for calculation of the top
#'   features.
#'
#' #transform and take return untransformed, top 5 features, and top 10 features
#' y <- transform(cc, dimReduce="var", nVarDims=c(NA, 5, 10))
#' names(y)
#'
# old .transData function for the variance filter part.
#
#   ###################
#   ###Dim Reduction
#   ###################
#   ##Check user inputs
#   ###################
#   #check valid options for dimReduce
#   varValues <- c("var","mad","cv")
#   if(any(!dimReduce %in% c("none","PCA",varValues)))
#     stop("invalid options for 'dimReduce' must be one of: 'none','PCA',",paste(varValues,collapse=","))
#
#   if(any(dimReduce!="none")){
#
#     ##Function to check and interpret values given
#     checkValues <- function(name){
#       ndims <- switch(as.character(name %in% varValues),
#                       "TRUE"=nVarDims, "FALSE"=nPCADims)
#       red <- dimReduce
#       if(any(is.na(ndims)) & name %in% red){ #if NA in ndims
#         if(length(ndims)==1){ #ndims is only a NA value -- assume user goofed and meant to do just "none"
#           if(length(red)==1) red<-"none"
#           if(length(red)>1) red <- red[-match(name, red)]
#         }
#         else{# otherwise user meant to do none *as well* as dimReduce with other values.
#           red <- unique(c("none", red)) #add 'none' and remove NA
#           ndims <- ndims[!is.na(ndims)]
#         }
#       }
#       dimReduce <<- red
#       if( name %in% varValues) nVarDims<<- ndims
#       if(name =="PCA") nPCADims <<- ndims
#     }
#
#
#     lapply(c("PCA", varValues), checkValues)
#
#     dimReduce <- unique(dimReduce)
#     nVarDims <- unique(nVarDims)
#     nPCADims <- unique(nPCADims)
#
#     xPCA <- xVAR <- xNone <-NULL #possible values
#     #logical as to whether return single matrix or list of matrices
#     listReturn<- !(length(dimReduce)==1 &&
#                      (dimReduce=="none" ||
#                         (dimReduce=="PCA" & length(nPCADims)==1) ||
#                         (dimReduce %in% varValues & length(nVarDims)==1)))
#     whFeatures <- NULL
#
#     xCL <- x
#     if(!is.null(clustering)){
#       if(any(!is.numeric(clustering)))
#         stop("clustering vector must be numeric")
#       if(length(clustering)!=ncol(x))
#         stop("clustering must be vector of length equal to columns of x")
#       if(all(clustering<0))
#         stop("All entries of clustering are negative")
#       if(sum(clustering<0)==ncol(x)-1)
#         stop("only one value in clustering not negative, cannot do dim reduction")
#       if(any(clustering<0))
#         xCL<-x[, -which(clustering<0)]
#     }
#
#
#     ##################
#     #PCA dim reduction
#     ##################
#
#     if("PCA" %in% dimReduce){
#
#       ######Check dimensions
#       if(max(nPCADims)>NROW(x))
#         stop("the number of PCA dimensions must be strictly less than the number of rows of input data matrix")
#       if(min(nPCADims)<=0)
#         stop("the number of PCA dimensions must be a value greater than 0")
#
#       pctReturn <- any(nPCADims < 1)
#       if(max(nPCADims)>100)
#         warning("the number PCA dimensions to be selected is greater than 100. Are you sure you meant to choose to use PCA dimensionality reduction rather than the top most variable features?")
#
#       ######Check zero variance genes:
#       rowvars <- matrixStats::rowVars(x)
#       if(any(rowvars==0)) {
#         if(all(rowvars==0)) {
#           stop("All features have zero variance.")
#         }
#         warning("Found features with zero variance.\nMost likely these are features with 0 across all samples.\nThey will be removed from PCA dimensionality reduction step.")
#       }
#       if(pctReturn) {
#         prcObj<-stats::prcomp(t(x[which(rowvars>0),]),center=TRUE,scale=TRUE)
#         prvar<-prcObj$sdev^2 #variance of each component
#         prvar<-prvar/sum(prvar)
#         prc<-t(prcObj$x)
#       } else {
#         prc <- t(.pca(t(x[which(rowvars>0),]), center=TRUE, scale=TRUE,
#                     k=max(nPCADims)))
#       }
#
#       if(pctReturn & NCOL(prc) != NCOL(origX))
#         stop("error in coding of principal components.")
#
#       if(any(nPCADims > NCOL(prc)))
#         stop("error in coding of principal components.")
#
#       if(!listReturn){ #nPCADims length 1; just return single matrix
#         if(pctReturn) {
#           nPCADims <- which(cumsum(prvar)>nPCADims)[1] #pick first pca coordinate with variance > value
#           xRet <- prc[seq_len(nPCADims),]
#         } else {
#           xRet <- prc
#         }
#       } else{
#         if(pctReturn){
#           whPct <- which(nPCADims<1)
#           pctNDims <- sapply(nPCADims[whPct], function(pct){
#             val<-which(cumsum(prvar)>pct)[1] #pick first pca coordinate with variance > value
#             if(length(val)==0) val<-length(prvar) #in case some numerical problem
#             return(val)
#           })
#           #if(any(is.na(pctNDims))) browser()
#           nPCADims[whPct]<-pctNDims
#         }
#         xPCA <- lapply(nPCADims,function(nn){prc[seq_len(nn),]})
#         names(xPCA)<-paste("nPCAFeatures=",nPCADims,sep="")
#       }
#     }
#
#     ##################
#     #Feature variability dim reduction
#     ##################
#     #for each dim reduction method requested
#     capwords <- function(s, strict = FALSE) { #From help of tolower
#       cap <- function(s) paste(toupper(substring(s, 1, 1)),
#                                {s <- substring(s, 2); if(strict) tolower(s) else s},
#                                sep = "", collapse = " " )
#       sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
#     }
#     doVarReduce<-function(name){
#       fun<-switch(name,"var"=stats::var,"mad"=stats::mad,"cv"=function(x){stats::sd(x)/mean(x)})
#
#       if(name %in% dimReduce){
#         if(max(nVarDims)>NROW(xCL)) stop("the number of most variable features must be strictly less than the number of rows of input data matrix")
#         if(min(nVarDims)<1) stop("the number of most variable features must be equal to 1 or greater")
#         if(min(nVarDims)<50 & NROW(xCL)>1000) warning("the number of most variable features to be selected is less than 50. Are you sure you meant to choose to use the top most variable features rather than PCA dimensionality reduction?")
#         varX<-apply(xCL,1,fun)
#         ord<-order(varX,decreasing=TRUE)
#         xVarOrdered<-x[ord,]
#         if(NCOL(xVarOrdered)!=NCOL(origX)) stop("error in coding of most variable.")
#         if(!listReturn){ #just return single matrix
#           xRet<-xVarOrdered[1:nVarDims,]
#           whFeatures<-ord[1:nVarDims]
#           return(list(x=xRet,whFeatures=whFeatures))
#         }
#         else{ #otherwise make it a list
#           xLIST<-lapply(nVarDims,function(nn){xVarOrdered[1:nn,]})
#           listName<-paste("n",toupper(name),"Features=",sep="")
#           names(xLIST)<-paste(listName,nVarDims,sep="")
#           return(xLIST)
#         }
#       }
#       else return(NULL)
#     }
#     if(any(dimReduce %in% varValues)){
#       dimReduceVar<-dimReduce[dimReduce %in% varValues]
#       # browser()
#       if(!listReturn & length(dimReduceVar)==1){
#         out<-doVarReduce(dimReduce)
#         xRet<-out$x
#         whFeatures<-out$whFeatures
#       }
#       else{
#         varOut<-lapply(dimReduceVar,doVarReduce)
#         xVAR<-unlist(varOut,recursive=FALSE)
#       }
#     }
#     #     if("var" %in% dimReduce & all(!is.na(nVarDims))){ #do PCA dim reduction
#     #       if(max(nVarDims)>NROW(x)) stop("the number of most variable features must be strictly less than the number of rows of input data matrix")
#     #       if(min(nVarDims)<1) stop("the number of most variable features must be equal to 1 or greater")
#     #       if(min(nVarDims)<50 & NROW(x)>1000) warning("the number of most variable features to be selected is less than 50. Are you sure you meant to choose to use the top most variable features rather than PCA dimensionality reduction?")
#     #       varX<-apply(x,1,mad)
#     #       ord<-order(varX,decreasing=TRUE)
#     #       xVarOrdered<-x[ord,]
#     #       if(NCOL(xVarOrdered)!=NCOL(origX)) stop("error in coding of principle components.")
#     #       if(length(nVarDims)==1 & length(dimReduce)==1){ #just return single matrix
#     #         x<-xVarOrdered[1:nVarDims,]
#     #         whFeatures<-ord[1:nVarDims]
#     #
#     #       }
#     #       else{ #otherwise make it a list
#     #         xVAR<-lapply(nVarDims,function(nn){xVarOrdered[1:nn,]})
#     #         names(xVAR)<-paste("nVarFeatures=",nVarDims,sep="")
#     #         listReturn<-TRUE
#     #       }
#     #     }
#     if("none" %in% dimReduce){
#       if(listReturn) xNone<-list("noDimReduce"=x)
#       else xRet<-x
#     }
#
#   }
#   else{
#     listReturn<-FALSE
#     whFeatures<-NULL
#     xRet<-x
#
#   }
#
#
#   if(listReturn) xRet<-c(xNone,xVAR,xPCA)
#   return(list(x=xRet,transFun=transFun,whMostVar=whFeatures))
# }