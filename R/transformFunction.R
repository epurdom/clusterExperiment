#' Transform the original data in a ClusterExperiment object
#'
#' Provides the transformed data (as defined by the object), as well as
#' dimensionality reduction.
#'
#' @param x a ClusterExperiment object.
#' @param nPCADims Numeric vector giving the number of PC dimensions to use in
#'   PCA dimensionality reduction. If NA no PCA dimensionality reduction is
#'   done.
#' @param nVarDims Numeric vector giving the number of features (e.g. genes) to
#'   keep, based on MAD variability.
#' @param dimReduce Character vector specifying the dimensionality reduction to
#'   perform, any combination of 'none', 'PCA' and 'mostVar'. See details.
#'
#' @details The data matrix defined by \code{assay(x)} is transformed based on
#'   the transformation function defined in x. If \code{dimReduce="none"} the
#'   transformed matrix is returned. Otherwise, the user can request
#'   dimensionality reduction of the transformed data via \code{dimReduce}.
#'   'PCA' refers to PCA of the transformed data with the top nPCADims kept.
#'   'mostVar' refers to keeping the top most variable features (defined by
#'   taking the MAD across all samples), and nVarDims defines how many such
#'   features to keep.
#' @details The PCA uses prcomp on \code{t(assay(x))} with \code{center=TRUE}
#'   and \code{scale=TRUE} (i.e. the feature are centered and scaled), so that
#'   it is performing PCA on the correlation matrix of the features.
#'   \code{dimReduce}, \code{nPCADims}, \code{nVarDims} can all be a vector of
#'   values, in which case a list will be returned with the appropriate datasets
#'   as elements of the list.
#'
#' @return If \code{dimReduce}, \code{nPCADims}, \code{nVarDims} are all of
#'   length 1, a matrix will be returned of the same dimensions as
#'   \code{assay(x)}. If these terms are vectors, then a list of data matrices
#'   will be return, each corresponding to the multiple choices implied by these
#'   parameters.
#'
#' @importFrom matrixStats rowVars
#'
#' @examples
#' mat <- matrix(data=rnorm(200), ncol=10)
#' mat[1,1] <- -1 #force a negative value
#' labels <- gl(5, 2)
#'
#' cc <- clusterExperiment(mat, as.numeric(labels), transformation =
#' function(x){x^2}) #define transformation as x^2
#'
#' #transform and take top 3 dimensions
#' x <- transform(cc, dimReduce="PCA", nPCADims=3)
#'
#' #transform and take return untransformed, top 5 features, and top 10 features
#' y <- transform(cc, dimReduce="mostVar", nVarDims=c(NA, 5, 10))
#' names(y)
#'
#' z<-transform(cc) #just return tranformed data
#' @export
#' @name transform
#' @aliases transform,ClusterExperiment-method
setMethod(
  f = "transform",
  signature = "ClusterExperiment",
  definition = function(x,nPCADims=NA,nVarDims=NA,dimReduce="none") {
    fun<-transformation(x)
    dat<-assay(x)
    return(.transData(dat,transFun=fun,nPCADims=nPCADims,nVarDims=nVarDims,dimReduce=dimReduce)$x)
  }
)

#function to transform assay data into clustering data (or other normal-like data input)
#Note for developers:
# .transData (unlike transform() ) returns a list:
# 1st element is the transformed data
# if npcs=NA or length of npcs=1, transformed data is matrix; otherwise returns list of data matrices.
# 2nd element is the transformation function
# The 2nd element is useful if function allows user to say isCount=TRUE so you can then actually get the transformation function out for defining ClusterExperiment Object)
# 3rd element is the index of most variable (if dimReduce="mostVar" and returns a simple matrix) otherwise NULL
.transData<-function(x,transFun=NULL,isCount=FALSE,nPCADims,nVarDims,dimReduce)
{
  origX<-x
  #transform data
  if(is.null(transFun)){
    transFun<-if(isCount) function(x){log2(x+1)} else function(x){x}
  }
  x<-try(transFun(x),silent=TRUE)
  if(inherits(x, "try-error")) stop(paste("User-supplied `transFun` produces error on the input data matrix:\n",x))
  if(any(is.na(x))) stop("User-supplied `transFun` produces NA values")
  #browser()
  ###################
  ###Dim Reduction
  ###################
  ##Check user inputs
  ###################
  listReturn<-FALSE
  whFeatures<-NULL
  #check valid options for dimReduce
  if(any(!dimReduce %in% c("none","PCA","mostVar"))) stop("invalid options for 'dimReduce' must be one of: 'none','PCA',or 'mostVar'")
  if(any(dimReduce!="none")){
    if(any(is.na(nPCADims)) & "PCA" %in% dimReduce){
      if(length(nPCADims)==1){
        if(length(dimReduce)==1) dimReduce<-"none" #assume user goofed and meant to do none
        if(length(dimReduce)>1) dimReduce<-dimReduce[-match("PCA",dimReduce)] #assume user goofed and didn't mean to also include
      }
      else{
        #add 'none' option to dimReduce and get rid of NA
        dimReduce<-unique(c("none",dimReduce))
        nPCADims<-nPCADims[!is.na(nPCADims)]
      }
    }
    if(any(is.na(nVarDims)) & "mostVar" %in% dimReduce){
      if(length(nVarDims)==1){
        if(length(dimReduce)==1) dimReduce<-"none" #assume user goofed and meant to do none
        if(length(dimReduce)>1) dimReduce<-dimReduce[-match("mostVar",dimReduce)] #assume user goofed and didn't mean to also include

      }
      else{#assume user meant to do none as well as dimReduce with other values.
        dimReduce<-unique(c("none",dimReduce)) #add 'none' and remove NA
        nVarDims<-nVarDims[!is.na(nVarDims)]
      }
    }
    dimReduce<-unique(dimReduce)
    nVarDims<-unique(nVarDims)
    nPCADims<-unique(nPCADims)
    xPCA<-xVAR<-xNone<-NULL #possible values
    #for each dim reduction method requested
    if("PCA" %in% dimReduce & !all(is.na(nPCADims))){ #do PCA dim reduction
      if(max(nPCADims)>NROW(x)) stop("the number of PCA dimensions must be strictly less than the number of rows of input data matrix")
      if(min(nPCADims)<1) stop("the number of PCA dimensions must be equal to 1 or greater")
      if(max(nPCADims)>100) warning("the number PCA dimensions to be selected is greater than 100. Are you sure you meant to choose to use PCA dimensionality reduction rather than the top most variable features?")

      rowvars <- matrixStats::rowVars(x)
      if(any(rowvars==0)) {

        if(all(rowvars==0)) {
          stop("All features have zero variance.")
        }

        warning("Found features with zero variance.\nMost likely these are features with 0 across all samples.\nThey will be removed from PCA.")
      }

      prc<-t(stats::prcomp(t(x[rowvars>0,]),center=TRUE,scale=TRUE)$x)
      if(NCOL(prc)!=NCOL(origX)) stop("error in coding of principle components.")
      if(length(nPCADims)==1 & length(dimReduce)==1){ #just return single matrix
        x<-prc[1:nPCADims,]

      }
      else{
        xPCA<-lapply(nPCADims,function(nn){prc[1:nn,]})
        names(xPCA)<-paste("nPCAFeatures=",nPCADims,sep="")
        listReturn<-TRUE
      }
    }
    if("mostVar" %in% dimReduce & all(!is.na(nVarDims))){ #do PCA dim reduction
      if(max(nVarDims)>NROW(x)) stop("the number of most variable features must be strictly less than the number of rows of input data matrix")
      if(min(nVarDims)<1) stop("the number of most variable features must be equal to 1 or greater")
      if(min(nVarDims)<50 & NROW(x)>1000) warning("the number of most variable features to be selected is less than 50. Are you sure you meant to choose to use the top most variable features rather than PCA dimensionality reduction?")
      varX<-apply(x,1,mad)
      ord<-order(varX,decreasing=TRUE)
      xVarOrdered<-x[ord,]
      if(NCOL(xVarOrdered)!=NCOL(origX)) stop("error in coding of principle components.")
      if(length(nVarDims)==1 & length(dimReduce)==1){ #just return single matrix
        x<-xVarOrdered[1:nVarDims,]
        whFeatures<-ord[1:nVarDims]

      }
      else{ #otherwise make it a list
        xVAR<-lapply(nVarDims,function(nn){xVarOrdered[1:nn,]})
        names(xVAR)<-paste("nVarFeatures=",nVarDims,sep="")
        listReturn<-TRUE
      }
    }
    if("none" %in% dimReduce & length(dimReduce)>1){
      xNone<-list("noDimReduce"=x)
      listReturn<-TRUE
    }

  }
  #browser()

  if(listReturn) x<-c(xNone,xVAR,xPCA)
  return(list(x=x,transFun=transFun,whMostVar=whFeatures))
}
