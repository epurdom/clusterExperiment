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
	labs<-clusterLabels(x)[wh]
	listByClusters<-strsplit(labs,",")
	
	listByParams<-sapply(1:length(listByClusters[[1]]),function(kk){
		sapply(listByClusters,.subset2,kk)})
	nameList<-lapply(1:ncol(listByParams),function(kk){
		yy<-listByParams[,kk]
		n<-sapply(strsplit(yy,"="),.subset2,1)
	})
	
	
	vals<-lapply(1:ncol(listByParams),function(kk){
		yy<-listByParams[,kk]
		v<-sapply(strsplit(yy,"="),function(zz){if(length(zz)>1).subset2(zz,2) else zz}) #deal with problem if no "=", e.g. 'noDimReduce' ...
		numV<-suppressWarnings(as.numeric(v))
		if(any(is.na(numV))) return(v) else return(numV)
	})
	
	##Deal with dim Reduce that can be multiple labels
	##This has to be manually updated when have new values in transform!
	dimChoices<-toupper(c("PCA","var","mad","cv"))
	dimValues <- c("noDimReduce",paste("n",dimChoices,"Features",sep=""))
	whDimReduce<-which(sapply(nameList,function(yy){any(yy %in% dimValues)}))
	if(length(whDimReduce)>1) stop("coding error: not expecting to have more than one dimReduce parameter")
	if(length(whDimReduce)==1){
		if(length(unique(nameList[[whDimReduce]]))==1){ 
			#only single value for dimReduce, so works like normal; use standard code
			nms<-sapply(nameList,unique)
			if(nms[[whDimReduce]]=="noDimReduce"){
				#if no Dim reduction on all, then just remove (would this even show up in clusterLabel)
				nms<-nms[-whDimReduce]
				vals<-vals[-whDimReduce]	
			}
		}
		else{
			##First, create variable giving type of dimReduce
			dimReduceType<-nameList[[whDimReduce]]
			nameList[[whDimReduce]]<-"nDimReduce"
			nms<-sapply(nameList,function(n){
				n<-unique(n)
				if(length(n)>1) stop(paste("found the following multiple values for the parameter name:",paste(n,collapse=",")))
				return(n)
			})
			##Fix up any that had 'noDimReduce' to be NA in vals and then make it numeric
			if(any(dimReduceType=="noDimReduce")){
				whNoDim<-which(dimReduceType=="noDimReduce")
				vals[[whDimReduce]][whNoDim]<-NA
				vals[[whDimReduce]]<-as.numeric(vals[[whDimReduce]])
			}
			##add dimReduceType to values
			dimReduceType<-gsub("Features","",dimReduceType)
			for(dimMethod in dimChoices){
				dimReduceType[dimReduceType==paste("n",dimMethod,sep="")]<-dimMethod
			}
			dimReduceType[dimReduceType=="noDimReduce"]<-"none"
			
			if(whDimReduce==1){
				nms<-c("dimReduce",nms)
				vals<-c(list(dimReduceType),vals)
			}
			else{
				nms<-c(nms[1:(whDimReduce)],"typeDimReduce",nms[whDimReduce:length(nms)])
				vals<-c(vals[1:(whDimReduce)],list(dimReduceType),vals[whDimReduce:length(nms)])
			}
		}
	
	}
	else{
		nms<-sapply(nameList,function(n){
			n<-unique(n)
			if(length(n)>1) stop(paste("found the following multiple values for the parameter name:",paste(n,collapse=",")))
			return(n)
		})
	}
	
	
	
	vals<-do.call("data.frame",c(vals,list(stringsAsFactors=FALSE)))
	colnames(vals)<-nms
	vals<-data.frame("clusteringIndex"=wh,vals)
	return(vals)
}
)