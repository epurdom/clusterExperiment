#' @title Return matrix from ClusterExperiment with reduced dimensions 
#' @description Returns a matrix of data from a \code{ClusterExperiment} object
#'   based on the choices of dimensionality reduction given by the user.
#' @rdname reduceFunctions
#' @aliases getReducedData
#' @param object  For \code{makeReducedDims},\code{makeFilterStats}, \code{defaultNDims} either matrix-like, \code{SingleCellExperiment}, or \code{ClusterExperiment} object. For \code{getReducedData} only a \code{ClusterExperiment} object allowed.
#' @param nDims The number of dimensions to keep from \code{reduceMethod}. If
#'   missing calls \code{\link{defaultNDims}}.
#' @param whichCluster an integer index or character string that identifies
#'   which cluster should be used to make the dendrogram. Default is
#'   primaryCluster.
#' @param whichAssay numeric or character specifying which assay to use. See
#'   \code{\link[SummarizedExperiment]{assay}} for details.
#' @param filterIgnoresUnassigned logical. Whether filtering statistics should 
#'   ignore the unassigned samples within the clustering. Only relevant if 
#'   'reduceMethod' matches one of built-in filtering statistics in 
#'   \code{\link{listBuiltInFilterStats}()}, in which case the clustering 
#'   identified in \code{whichCluster} is passed to \code{makeFilterStats} and
#'   the unassigned samples are excluded in calculating the statistic. See 
#'   \code{\link{makeFilterStats}}  for more details. 
#' @param returnValue The format of output. Users will generally want to keep
#'   the default (see details)
#' @param reducedDimName The name given to the reducedDims slot storing result 
#'   (if \code{returnValue="object"}). If missing, the function will create a
#'   default name: if \code{reduceMethod} is a dimensionality reduction, then
#'   \code{reduceMethod} will be given as the name; if a filtering statistic,
#'   "filteredBy_" followed by \code{reduceMethod}.
#' @details This function determines the matrix of values that can be used for
#'   computation based on the user's choice of dimensionality methods. The
#'   methods can be either of the filtering kind or the more general
#'   dimensionality reduction. The function will first look at any stored
#'   \code{ReducedDims} or filtering statistics already present in the data, and
#'   if missing, will assume that \code{reduceMethod} is one of the built-in
#'   method provided by the package and calculate the necessary. Note that if
#'   \code{reduceMethod} is a filtering statistic, in addition to filtering the
#'   features, the function will also perform the stored transformation of the
#'   data.
#' @details Note that this is used internally by functions, but is mainly only
#'   of interest for the user if they want to have the filtered, transformed
#'   data available as a matrix for continual use.
#' @details If \code{returnValue="object"}, then the output is a single, updated
#'   \code{ClusterExperiment} object with the reduced data matrix stored as an 
#'   element of the list in \code{reducedDims} slot (with name given by
#'   \code{reducedDimName} if given). If "list", then a list with one element
#'   that is the object and the other that is the reduced data matrix. Either
#'   way, the object returned in the list will be updated to contain with the
#'   filtering statistics or the dimensionality reduction. The only difference
#'   is that if "list", the reduced dimension matrix is NOT saved in the object
#'   (and so only really makes a difference if the \code{reduceMethod} argument
#'   is a filtering method). The option "list" is mainly for internal use, where
#'   we do not want to continually save subseted datasets.
#' @details If \code{nDims} is missing, it will be given a default value 
#'   depending on the value of \code{reduceMethod}. See 
#'   \code{\link{defaultNDims}} for details.
#' @details If \code{filterIgnoresUnassigned} is missing, then it is set to TRUE
#'   \emph{unless}: \code{reduceMethod} matches a stored filtering statistic in
#'   \code{rowData} AND does not match a built-in filtering method provided by
#'   the package.
#' @return If \code{returnValue="object"}, a \code{ClusterExperiment} object.
#' @return If \code{returnValue="list"} a list with elements:
#' \itemize{
#' \item{\code{objectUpdate}}{object, potentially updated if had to calculate dimensionality reduction or filtering statistic}
#' \item{\code{dataMatrix}}{the reduced dimensional matrix with the samples in columns, features in rows}
#' }
#' @seealso \code{\link{makeFilterStats}},\code{\link{makeReducedDims}}, \code{\link{filterData}}, \code{\link[SingleCellExperiment]{reducedDim}}
setMethod(
  f = "getReducedData",
  signature = "ClusterExperiment",
  definition = function(object,reduceMethod,filterIgnoresUnassigned, 
                        nDims=defaultNDims(object,reduceMethod),whichCluster="primary", 
                        whichAssay=1, returnValue=c("object","list"),reducedDimName){
    if(isReducedDims(object,reduceMethod) & isFilterStats(object,reduceMethod)) stop(paste(reduceMethod,"is the name of both a stored filtering statistic and a stored dimensionality reduction -- cannot create reduced data"))
    whCl<-.convertSingleWhichCluster(object,whichCluster)
    returnValue<-match.arg(returnValue)
    reduceMethodName<-reduceMethod
    if(missing(filterIgnoresUnassigned)){
      if(isFilterStats(object,reduceMethod) & !isBuiltInFilterStats(reduceMethod)) filterIgnoresUnassigned<-FALSE
      else filterIgnoresUnassigned<-TRUE
    }
    if(filterIgnoresUnassigned & isBuiltInFilterStats(reduceMethod)){
      reduceMethodName<-.makeClusterFilterStats(reduceMethod,clusterLabels(object)[whCl])
    }
    if(length(nDims) > 1) {
      stop("getReducedData only handles one choice of dimensions.")
    }
    if(!is.na(nDims) & reduceMethod=="none") {
      warning("specifying nDims has no effect if reduceMethod==`none`")
    }
    
    #deal with name of new assay
    if(returnValue=="object" ){
      addToObject<-TRUE
      if(isReducedDims(object,reduceMethod)){
        addToObject<-FALSE 
        warning("will not add reduced dataset to object because already exists method with that name")
      }  
      else{
        if(missing(reducedDimName)){
          if(isFilterStats(object,reduceMethod) || isBuiltInFilterStats(reduceMethod)){
            reducedDimName<-paste0("filteredBy_",reduceMethodName)
          }
          else reducedDimName<-reduceMethod
        }
      }
    }
    else addToObject<-FALSE

    #check nDims less than nDims of saved object...otherwise will have to recalculate it, if possible, and give it new name
    if(isReducedDims(object,reduceMethod) && nDims > ncolReducedDims(object)[[reduceMethod]]){
        warning("requesting an existing dimensionality reduction but with greater number of dimensions than available. Will ignore nDims and use max dimensions")
        nDims<-ncolReducedDims(object)[reduceMethod]
    }
    ###Calculate filters/reduceMethod if needed...
    if(!isReducedDims(object,reduceMethod) & isBuiltInReducedDims(reduceMethod)){
      object<-makeReducedDims(object,reducedDims=reduceMethod,maxDims=nDims,whichAssay=whichAssay)
    }
    else if(!isFilterStats(object,reduceMethodName) & isBuiltInFilterStats(reduceMethod)){
      object<-makeFilterStats(object,filterStat=reduceMethod, whichAssay=whichAssay,
                              whichClusterIgnoreUnassigned=if(filterIgnoresUnassigned) whCl else NULL)
    }
    
    if(reduceMethod=="none")
      dat<-transformData(object, whichAssay=whichAssay)
    else if(isReducedDims(object,reduceMethod))
      dat<-t(reducedDim(object,type=reduceMethod)[,seq_len(nDims)])
    else if(isFilterStats(object,reduceMethodName))
      dat<-transformData(filterData(object,filterStats=reduceMethodName,percentile=nDims), whichAssay=whichAssay)
    else stop("'object' does not contain the given 'reduceMethod' value nor does 'reduceMethod' value match any built-in filters or dimensionality reduction options.")
    
    if(addToObject){
      if(isReducedDims(object,reduceMethod) & reduceMethod!=reducedDimName){
        #means was already calcualted and saved with name 'reduceMethod' but not right name
        wh<-which(reducedDimNames(object)==reduceMethod)
        if(length(wh)>1) stop("internal coding error -- multiple matches of reduceMethod")
        names(reducedDims(object))[wh]<-reducedDimName
      }
      else reducedDim(object,reducedDimName) <- t(dat)
    }
    if(returnValue=="list"){
      return(list(objectUpdate=object,dataMatrix=dat))
    }
    if(returnValue=="object"){
      return(object)
    }
  }
)
#redo filtering statistic on subset determined by cluster unassigned
.makeClusterFilterStats<-function(filterStats,clusterName){
  make.names(paste(filterStats,clusterName,sep="_"))
}
#reduce reduced dims to get larger number of Dims
.makeClusterReducedRedo<-function(reducedDims,nDims){
  make.names(paste(reducedDims,nDims,sep="_"))
}
#' @rdname reduceFunctions
#' @param reduceMethod character. A method (or methods) for reducing the size of
#'   the data, either by filtering the rows (genes) or by a dimensionality
#'   reduction method. Must either be 1) must match the name of a built-in
#'   method, in which case if it is not already existing in the object will be
#'   passed to \code{\link{makeFilterStats}} or \code{link{makeReducedDims}}, or
#'   2) must match a stored filtering statistic or dimensionality reduction in
#'   the object
#' @param typeToShow character (optional). If given, should be one of
#'   "filterStats" or "reducedDims" to indicate of the values in the
#'   reduceMethod vector, only show those corresponding to "filterStats" or
#'   "reducedDims" options.
#' @return \code{defaultNDims} returns a numeric vector giving the default
#'   dimensions the methods in \code{clusterExperiment} will use for reducing
#'   the size of the data. If \code{typeToShow} is missing, the resulting vector
#'   will be equal to the length of \code{reduceMethod}. Otherwise, it will be a
#'   vector with all the unique valid default values for the \code{typeToShow}
#'   (note that different dimensionality reduction methods can have different
#'   maximal dimensions, so the result may not be of length one in this case).
#' @details For a \code{reduceMethod} that corresponds to a filtering statistics
#'   the current default is 1000 (or the length of the number of features, if
#'   less). For a dimensionality reduction saved in the reducedDims slot the
#'   default is 50 or the maximum number of dimensions if less than 50.
#' @details \code{reduceMethod} will first be checked to see if it corresponds
#'   with an existing saved filtering statistic or a dimensionality reduction to
#'   determine which of these two types it is. If it does not match either, then
#'   it will be checked against the built in functions provided by the package.
#'   @examples 
#'   se<-SingleCellExperiment(matrix(rnorm(5000*100),nrow=5000,ncol=100))
#'   defaultNDims(se,"PCA")
#'   defaultNDims(se,"mad")
#' @aliases defaultNDims defaultNDims,SingleCellExperiment-method
#' 
setMethod( 
  f="defaultNDims",
  "SingleCellExperiment",
  function(object,reduceMethod,typeToShow){
    nDims<-rep(NA,length(reduceMethod))
		isFilter<-isBuiltInFilterStats(reduceMethod) | isFilterStats(object,reduceMethod)
		isRed<-isReducedDims(object,reduceMethod ) | isBuiltInReducedDims(reduceMethod)
    isAnyFilter<-any(isFilter)
    isAnyRed<-any(isRed)
    if(isAnyFilter)
      nDims[isBuiltInFilterStats(reduceMethod) | isFilterStats(object,reduceMethod)]<-min(1000,NROW(object))
    if(isAnyRed){
			if(any(isReducedDims(object,reduceMethod )))
      nDims[isReducedDims(object,reduceMethod)] <- ncolReducedDims(object)[reduceMethod[isReducedDims(object,reduceMethod )]]
    	if(any(!isReducedDims(object,reduceMethod )& isBuiltInReducedDims(reduceMethod))) nDims[!isReducedDims(object,reduceMethod )& isBuiltInReducedDims(reduceMethod)]<-min(c(50,dim(object)))
		}
		if(!missing(typeToShow)){ #means pick a single one for each type
      if(typeToShow=="filterStats"){
        if(any(isFilter)) nDims<-min(unique(nDims[isFilter])) else nDims<-NA
      }
      if(typeToShow=="reducedDims"){
        if(any(isRed)) nDims<-min(unique(nDims[isRed])) else nDims<-NA
      }
      if(length(nDims)==0) nDims<-NA
    }
    return(nDims)
    
  })


setMethod( 
  f="defaultNDims","matrixOrHDF5",function(object,...){
    return(defaultNDims(SingleCellExperiment(object),...))
  })

