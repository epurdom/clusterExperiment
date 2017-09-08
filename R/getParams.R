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
#' @details The method simply parses the \code{clusterLabels} of the indicated
#'   clusterings, relying on the specific format used by \code{clusterMany} to
#'   create labels. The function will only allow the parsing to be performed on
#'   those clusterings with a 'clusterMany' clusterType. If the user has
#'   manipulated the clusterLabels manually or manually identified the
#'   clusterType of a clustering as 'clusterMany', this function may create
#'   unexpected results or errors.
#' @details Specifically, it splits the label of each clustering by the
#'   character ",", as indicating the different parameters; this should return a
#'   value of form "ABC=123". The function then pulls out the numeric value
#'   ('123') and associates that value as the value of the parameter ('ABC')
#' @return Returns a data.frame where the column names are the parameter names,
#'   and the entries are the values of the parameter for the indicated
#'   clustering. The column 'clusteringIndex' identifies the index of the
#'   clustering in the full set of clusterings of the given ClusterExperiment
#'   object.
#' @export
setMethod(f = "getClusterManyParams",
          signature = signature(x = "ClusterExperiment"),
    function(x,whichClusters="clusterMany") {
					
	allClTypes<-clusterTypes(x)
	if(!"clusterMany" %in% allClTypes) stop("x does not have any clusterings that have clusterType exactly equal to 'clusterMany' ")	 #note to self: this makes it impossible to use function on old clusterMany results that have type clusterMany.1 or what ever.		
	wh <- .TypeIntoIndices(x, whClusters=whichClusters)
	if(length(wh)==0 || !"clusterMany" %in% allClTypes[wh]){
		warning("argument whichClusters either did not return any clusters or did not return any that were from 'clusterMany' type; using all clusterMany clusters")
		wh<-which(allClTypes=="clusterMany")
	}
	else{
		if(any(allClTypes[wh] != "clusterMany")){
			warning("some clusters indicated in 'whichClusters' do not have type 'clusterMany' and will be ignored.") #note to self: again, this makes it impossible to use function on old clusterMany results that have type clusterMany.1 or what ever.
			wh<-wh[allClTypes[wh] == "clusterMany"]
		}
		
	}
	labs<-clusterLabels(x)[wh]
	listByClusters<-strsplit(labs,",")
	
	listByParams<-sapply(1:length(listByClusters[[1]]),function(kk){
		sapply(listByClusters,.subset2,kk)})
	nms<-apply(listByParams,2,function(yy){
		n<-unique(sapply(strsplit(yy,"="),.subset2,1))
		if(length(n)>1) stop(paste("found the following multiple values for the parameter name:",paste(n,collapse=",")))
	})
	
	vals<-lapply(1:ncol(listByParams),function(kk){
		yy<-listByParams[,kk]
		v<-(sapply(strsplit(yy,"="),.subset2,2))
		numV<-as.numeric(v)
		if(any(is.na(numV))) return(v) else numV
	})
	vals<-do.call("data.frame",vals)
	colnames(vals)<-nms
	vals<-data.frame("clusteringIndex"=wh,vals)
	return(vals)
}
)