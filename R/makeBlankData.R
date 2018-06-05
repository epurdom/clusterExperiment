#' @rdname plottingFunctions
#' @name plottingFunctions
#' @aliases makeBlankData
#' @param data matrix with samples on columns and features on rows.
#' @param groupsOfFeatures list, with each element of the list containing a
#'   vector of numeric indices of features (rows).
#' @param groupsOfSamples list, with each element of the list containing a
#'   vector of numeric indices of samples (columns).
#' @param nBlankFeatures the number of blank lines to add in the data matrix to
#'   separate the groups of feature indices (will govern the amount of white space if
#'   data is then fed to heatmap.)
#' @param nBlankSamples the number of blank lines to add in the data matrix to
#'   separate the groups of sample indices (will govern the amount of white space if
#'   data is then fed to heatmap.)
#' @details \code{makeBlankData} pulls the data corresponding to the row indices
#'   in \code{groupsOfFeatures} adds lines of NA values into data between these
#'   groups. When given to heatmap, will create white space between these groups
#'   of features.
#' @seealso \code{\link{plotHeatmap}}
#' @return \code{makeBlankData} returns a list with items
#'  \itemize{
#'  \item{"dataWBlanks"}{ The data with the rows of NAs separating the given
#'  indices.}
#'  \item{"rowNamesWBlanks"}{ A vector of characters giving the rownames for the
#'  data, including blanks for the NA rows. These are not given as rownames to
#'  the returned data because they are not necessarily unique. However, they can be 
#'  given to the \code{labRow} argument of \code{\link[NMF]{aheatmap}} or
#'  \code{\link{plotHeatmap}}.}
#'  \item{"groupNamesWBlanks"}{ A vector of characters of the same length 
#'	as the number of rows of the new data (i.e. with blanks) giving the group name  
#'  for the data, indicating which group (i.e. which element of \code{groupsOfFeatures}
#'   list) the feature came from. If \code{groupsOfFeatures} has unique names, these 
#'   names will be used, other wise "Group1", "Group2", etc. The NA rows are given 
#'   NA values. 
#'   }	
#' }
#'
#' @export
#'
#' @examples
#' data(simData)
#'
#' x <- makeBlankData(simData[,1:10], groupsOfFeatures=list(c(5, 2, 3), c(20,
#' 34, 25)))
#' plotHeatmap(x$dataWBlanks,clusterFeatures=FALSE)
makeBlankData <- function(data,groupsOfFeatures=NULL,groupsOfSamples=NULL,nBlankFeatures = 1, nBlankSamples=1) {
	###Checks on input arguments
	if(is.null(groupsOfFeatures)& is.null(groupsOfSamples))
		stop("One of 'groupsOfFeatures' and 'groupsOfSamples' must be non-NULL")
  if (!is.null(groupsOfFeatures) && !is.list(groupsOfFeatures))
    stop("groupsOfFeatures must be a list")
  if (!is.null(groupsOfSamples) && !is.list(groupsOfSamples))
    stop("groupsOfSamples must be a list")

  if(!is.null(groupsOfFeatures)){
		testFeatureIndices <- sapply(groupsOfFeatures,function(x) {
    	is.numeric(x) && all(x %in% seq_len(NROW(data)))
  	})
	  if (!all(testFeatureIndices))
	    stop("Invalid list of indices in groupsOfFeatures.")
	}
	
  if(!is.null(groupsOfSamples)){
		testSampleIndices <- sapply(groupsOfSamples,function(x) {
	    is.numeric(x) && all(x %in% seq_len(NROW(data)))
	  })
	  if (!all(testSampleIndices))
	    stop("Invalid list of indices in groupsOfFeatures.")
	}

  if(is.null(rownames(data))) row.names(data)<-as.character(seq_len(nrow(data)))
	if(is.null(colnames(data))) colnames(data)<-as.character(seq_len(ncol(data)))
	rnames<-rownames(data)
	cnames<-colnames(data)
	####
	#Features
	######
	if(!is.null(groupsOfFeatures)){
	  #make list of data of feature groups
	  dataList <- lapply(groupsOfFeatures,function(ii){data[ii,,drop=FALSE]})
	  rnames <- lapply(groupsOfFeatures,function(ii){rnames[ii]})
	  #remove rownames to avoid warning
	  dataList<-lapply(dataList,function(x){
		  row.names(x)<-NULL
		  return(x)
	  })
	  #add NA rows between groups
	  naData <- matrix(NA,nrow = nBlankFeatures,ncol = ncol(data))
	  dataListMinus <- lapply(dataList[-length(dataList)],function(x) {
		  rbind(x,naData)
	  })  
	  newData <-
	    data.frame(do.call("rbind",c(dataListMinus,dataList[length(dataList)])))

	  #make names for this
	  rnamesMinus <-
	    lapply(head(rnames,-1),function(x) {
	      c(x,rep("",nBlankFeatures))
	    })
	  rnames <- unname(c(unlist(rnamesMinus),rnames[[length(rnames)]]))
  
  
	  if(!is.null(names(groupsOfFeatures)) && length(unique(names(groupsOfFeatures)))==length(names(groupsOfFeatures)) ){
		  gNames<-names(groupsOfFeatures)
	  }
	  else gNames<-paste("Feature Group",seq_along(groupsOfFeatures),sep="")
	  groupNames<- lapply(seq_along(dataList),function(i){
		  gname<-rep(gNames[i],times=length(groupsOfFeatures[[i]]))
		  if(i!=length(dataList)) gname<-c(gname,rep(NA,nBlankFeatures))
		  return(gname)
	  })
	  featureGroupNames<-unname(unlist(groupNames))
		data<-newData
		
	}
	else{
		rnames<-NULL
		featureGroupNames<-NULL
	}
	####
	#Samples
	######
	if(!is.null(groupsOfSamples)){
	  #make list of data of feature groups
	  dataList <- lapply(groupsOfSamples,function(ii){data[,ii,drop=FALSE]})
	  cnames <- lapply(groupsOfSamples,function(ii){cnames[ii]})
	  #remove rownames to avoid warning
	  dataList<-lapply(dataList,function(x){
			colnames(x)<-NULL
					  return(x)
	  })
	  #add NA columns between groups
	  naData <- matrix(NA,ncol = nBlankSamples,nrow = nrow(data),dimnames=NULL)
	  dataListMinus <- lapply(dataList[-length(dataList)],function(x) {
		  cbind(x,naData)
	  })  
	  newData <- unname(data.frame(do.call("cbind",c(dataListMinus,dataList[length(dataList)]))))
	  #make names for this
	  cnamesMinus <-
	    lapply(head(cnames,-1),function(x) {
	      c(x,rep("",nBlankSamples))
	    })
	  cnames <- unname(c(unlist(cnamesMinus),cnames[[length(cnames)]]))
  
  
	  if(!is.null(names(groupsOfSamples)) && length(unique(names(groupsOfSamples)))==length(names(groupsOfSamples)) ){
		  gNames<-names(groupsOfSamples)
	  }
	  else gNames<-paste("Sample Group",seq_along(groupsOfSamples),sep="")
	  groupNames<- lapply(seq_along(dataList),function(i){
		  gname<-rep(gNames[i],times=length(groupsOfSamples[[i]]))
		  if(i!=length(dataList)) gname<-c(gname,rep(NA,nBlankSamples))
		  return(gname)
	  })
	  sampleGroupNames<-unname(unlist(groupNames))
		data<-newData
		
	}
	else{
		cnames<-NULL
		sampleGroupNames<-NULL
	}

  #can't set rownames because not unique values!
  #rownames(newData)<-rnames
  return(list(dataWBlanks = data,rowNamesWBlanks = rnames,colNamesWBlanks=cnames, featureGroupNamesWBlanks=featureGroupNames,sampleGroupNamesWBlanks=sampleGroupNames))

}
