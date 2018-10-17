#' @rdname getClusterManyParams
#' @title Get parameter values of clusterMany clusters
#'   
#' @description Takes an input a ClusterExperiment object and returns the
#'   parameter values used in creating the clusters that were created by
#'   'clusterMany'
#'   
#' @aliases getClusterManyParams
#'   
#' @param x a ClusterExperiment object that contains clusterings from running
#'   \code{\link{clusterMany}}.
#' @param searchAll logical, indicating whether all clusterings with 
#'   \code{clusterMany} label should be allowed (i.e. including those from
#'    previous ones with labels like \code{clusterMany.1}), or only limited 
#'    to those in most recent workflow (default).
#' @param simplify logical. Whether to simplify the output so as to remove 
#'   features that do not change across the clusterings.
#' @details The method simply parses the \code{clusterLabels} of the indicated
#'   clusterings, relying on the specific format used by \code{clusterMany} to
#'   create labels. The function will only allow the parsing to be performed on
#'   those clusterings with a 'clusterMany' clusterType. If the user has
#'   manipulated the clusterLabels manually or manually identified the
#'   clusterType of a clustering as 'clusterMany', this function may create
#'   unexpected results or errors. Similarly, it cannot be used on 'clusterMany'
#'   results from an old iteration (e.g. that have type 'clusterMany.1')
#' @details Specifically, it splits the label of each clustering by the
#'   character ",", as indicating the different parameters; this should return a
#'   value of form "ABC=123". The function then pulls out the numeric value
#'   ('123') and associates that value as the value of the parameter ('ABC')
#' @return Returns a data.frame where the column names are the parameter names,
#'   and the entries are the values of the parameter for the indicated
#'   clustering. The column 'clusteringIndex' identifies the index of the
#'   clustering in the full set of clusterings of the given ClusterExperiment
#'   object.
#' @inheritParams ClusterExperiment-methods
#' @examples
#' data(simData)
#'
#' cl <- clusterMany(simData, nReducedDims=c(5, 10, 50), reduceMethod="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(TRUE,FALSE),
#' removeSil=c(TRUE,FALSE))
#' getClusterManyParams(cl)
#' @export
setMethod(
  f = "getClusterManyParams",
  signature = signature(x = "ClusterExperiment"),
  function(x,whichClusters="clusterMany",searchAll=FALSE,simplify=TRUE) {
    
    allClTypes<-clusterTypes(x)
	if(searchAll) whCM<-grep("clusterMany",allClTypes)
	else whCM<-which(allClTypes=="clusterMany")
    if(length(whCM)==0) stop("x does not have any clusterings that have clusterType 'clusterMany' ")	
	wh <- getClusterIndex(x, whichClusters=whichClusters,noMatch="removeSilently")
    if(length(wh)==0){
      warning("argument whichClusters did not return any clusters; using all clusterMany clusters")
      wh<-whCM
    }
    if(!any(wh %in% whCM)){
      warning("argument whichClusters did not return any clusters of type 'clusterMany'; using all clusterMany clusters")
      wh<-whCM
    }
    else{
      if(any(!wh %in%  whCM)){
        warning("some clusters indicated in 'whichClusters' do not have type 'clusterMany' and will be ignored.") 
        wh<-wh[wh %in% whCM]
      }
      
    }
    infoList<-clusteringInfo(x)[wh]
    params<-do.call("rbind",lapply(infoList,.subset2,"choicesParam"))
    if(simplify){
      # notAllNA<-which(apply(params,2,function(x){!all(is.na(x))}))
      # params<-params[,notAllNA]
      
      notAllSame<-which(apply(params,2,function(z){length(unique(z))>1}))
      params<-params[,notAllSame,drop=FALSE]
    }
		row.names(params)<-clusterLabels(x)[wh]
    params<-data.frame("clusteringIndex"=wh,params)
    return(params)
  }
)

