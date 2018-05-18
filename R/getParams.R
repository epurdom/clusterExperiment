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
#' @param whichClusters The indices (or clusterLabels) of those clusters whose
#'   labels will be parsed to determine the parameters; should be subset of the
#'   clusterMany results.
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
  function(x,whichClusters="clusterMany",simplify=TRUE) {
    
    allClTypes<-clusterTypes(x)
    if(!"clusterMany" %in% allClTypes) stop("x does not have any clusterings that have clusterType exactly equal to 'clusterMany' ")	 #note to self: this makes it impossible to use function on old clusterMany results that have type clusterMany.1 or what ever.		
    wh <- if(is.character(whichClusters)).TypeIntoIndices(x, whClusters=whichClusters) else whichClusters
    if(length(wh)==0){
      warning("argument whichClusters did not return any clusters; using all clusterMany clusters")
      wh<-which(allClTypes=="clusterMany")
    }
    if(!"clusterMany" %in% allClTypes[wh]){
      warning("argument whichClusters did not return any clusters of type 'clusterMany'; using all clusterMany clusters")
      wh<-which(allClTypes=="clusterMany")
    }
    else{
      if(any(allClTypes[wh] != "clusterMany")){
        warning("some clusters indicated in 'whichClusters' do not have type 'clusterMany' and will be ignored.") #note to self: again, this makes it impossible to use function on old clusterMany results that have type clusterMany.1 or what ever.
        wh<-wh[allClTypes[wh] == "clusterMany"]
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

