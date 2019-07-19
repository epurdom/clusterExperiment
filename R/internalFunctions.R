.onUnload <- function (libpath) {
  library.dynam.unload("clusterExperiment", libpath)
}

.mynote<-function(x){
	message(paste("Note:",x))
}

# Small function to identify what type of CoClustering information is stored in co-clustering slot
# (might make sense to make it a slot, but then have to update object, for something not so important...)
.typeOfCoClustering<-function(ceObj){
	if(is.null(ceObj@coClustering)) return(NULL)
    if(is.null(dim(ceObj@coClustering))) return("indices")
        else if(!isSymmetric(ceObj@coClustering)) return("clusterings")
            else return("distance")
}

.addPrefixToClusterNames<-function(ceObj,prefix,whCluster){
  ceLegend<-clusterLegend(ceObj)[[whCluster]]
  cl<-ceLegend[,"clusterIds"]
  ceLegend[,"name"]<-numericalAsCharacter(values=cl,prefix=prefix)
  clusterLegend(ceObj)[[whCluster]]<-ceLegend
  return(ceObj)
}
#' @title Convert numeric values to character that sort correctly
#' @description Small function that takes as input integer values (or values
#'   that can be converted to integer values) and converts them into character
#'   values that are 'padded' with zeros at the beginning of the numbers so that
#'   they will sort correctly.
#' @param values vector of values to be converted into sortable character values
#' @param prefix optional character string that will be added as prefix to the
#'   result
#' @details The function determines the largest value and adds zeros to the
#'   front of smaller integers so that the resulting characters are the same
#'   number of digits. This allows standard sorting of the values to correctly
#'   sort.
#' @details The maximum number of zeros that will be added is 3. Input integers
#'   beyond that point will not be correctly fixed for sorting.
#' @details Negative integers will not be corrected, but left as-is
#' @return A character vector
#' @seealso \code{\link[stringr]{str_pad}}
#' @examples
#' numericalAsCharacter(c(-1, 5,10,20,100))
#' @importFrom stringr str_pad
#' @export
numericalAsCharacter<-function(values,prefix=""){
	if(is.factor(values)) values<-as.character(values)
	values<-suppressWarnings(as.integer(values))
	if(any(is.na(values))) stop("input must convert to numeric (integer) vector")
	whPos<-which(values >0)
	if(length(whPos)>0) largestLength<-max(values[whPos])
	values<-as.character(values)
	if(length(whPos)>0){
		pad<-if(largestLength<100) 2 else if(largestLength<1000) 3 else 4			
		values[whPos]<-paste(prefix,stringr::str_pad(values[whPos],width=pad,pad="0"),sep="")
  	} 
	return(values)
}

.addNewResult<-function(newObj,oldObj){
  retval<-addClusterings(newObj,oldObj) #want most recent addition on top of clusterMatrix
  #erases dendrogram from oldObj -- only keeps newObj -- so need to put it back if wasn't already there
  if(is.na(retval@dendro_index) & !is.na(newObj@dendro_index)) stop("Coding error -- addClusterings lost dendro_index")
  if(is.na(retval@merge_index) & !is.na(newObj@merge_index)) stop("Coding error -- addClusterings lost merge_index")
  if(is.na(retval@dendro_index) & !is.na(oldObj@dendro_index)){
    retval@dendro_samples<-oldObj@dendro_samples
    retval@dendro_clusters<-oldObj@dendro_clusters
    retval@dendro_index<-oldObj@dendro_index+nClusterings(newObj) #update index to where dendrogram from
  }
  if(is.na(retval@merge_index) & !is.na(oldObj@merge_index)){
    retval@merge_index<-oldObj@merge_index+nClusterings(newObj) #update index to where merge from
    retval@merge_nodeMerge<-oldObj@merge_nodeMerge
    retval@merge_cutoff<-oldObj@merge_cutoff
    retval@merge_method<-oldObj@merge_method
		retval@merge_demethod<-oldObj@merge_demethod
  }
  if(is.null(retval@merge_nodeProp) & !is.null(oldObj@merge_nodeProp)){
    retval@merge_nodeProp<-oldObj@merge_nodeProp
    retval@merge_dendrocluster_index<-oldObj@merge_dendrocluster_index+nClusterings(newObj) #update index to where merge from
  }
  #put back orderSamples, coClustering
  if(all(retval@orderSamples==seq_len(nSamples(retval))) & !all(oldObj@orderSamples==seq_len(nSamples(retval)))) retval@orderSamples<-oldObj@orderSamples
  if(is.null(retval@coClustering)) retval@coClustering<-oldObj@coClustering
  retval<-.addBackSEInfo(newObj=retval,oldObj=oldObj) #make sure keeps SE info
  #   Note: .addBackSEInfo calls ClusterExperiment (i.e. validates)
  return(retval)
}

#this function keeps everything from new, except grabs SE info from old
.addBackSEInfo<-function(newObj,oldObj){
  retval<-ClusterExperiment(as(oldObj,"SingleCellExperiment"),
                            clusters=clusterMatrix(newObj),
                            transformation=transformation(newObj),
                            clusterTypes=clusterTypes(newObj),
                            clusterInfo=clusteringInfo(newObj),
                            orderSamples=orderSamples(newObj),
                            coClustering=coClustering(newObj),
                            dendro_samples=newObj@dendro_samples,
                            dendro_clusters=newObj@dendro_clusters,
                            dendro_index=newObj@dendro_index,
                            merge_index=newObj@merge_index,
                            merge_cutoff=newObj@merge_cutoff,
                            merge_dendrocluster_index=newObj@merge_dendrocluster_index,
                            merge_nodeProp=newObj@merge_nodeProp,
                            merge_nodeMerge=newObj@merge_nodeMerge,
                            merge_method=newObj@merge_method,
                            merge_demethod=newObj@merge_demethod,
                            primaryIndex=primaryClusterIndex(newObj),
                            clusterLegend=clusterLegend(newObj),
                            checkTransformAndAssay=FALSE
  )
  return(retval)
}



#Returns NULL if no sample data
.pullColData<-function(ce,wh,fixNA=c("keepNA","unassigned","missing")){
  fixNA<-match.arg(fixNA)
  if(!is.null(wh)){
    sData<-colData(ce)
    if(!is.logical(wh)){
      
      if(NCOL(sData)==0) stop("no colData for object data, so cannot pull sample data")
      if(is.character(wh)){
        if(all(wh=="all")) wh<-seq_len(NCOL(sData))
        else{
          if(!all(wh %in% colnames(sData))) stop("Invalid names for pulling sample data (some do not match names of colData)")
          else wh<-match(wh,colnames(sData))
        }
      }
      else if(is.numeric(wh)){
        if(!all(wh %in% seq_len(NCOL(sData)))) stop("Invalid indices for for pulling sample data (some indices are not in 1:NCOL(colData)")
      }
      else stop("invalid values for pulling sample data from colData of object")
      sData<-as.data.frame(sData[,wh,drop=FALSE])
    }
    else{ #if 
      if(wh) sData<- colData(ce)
      else sData<-NULL
    }
  }
  else sData<-NULL
  if(!is.null(sData) && fixNA!="keepNA"){
    newValue<-switch(fixNA,"unassigned"=-1,"missing"=-2)
    fixNAFunction<-function(x,newValue){
      if(is.factor(x)){ #change to character
        waslevels<-levels(x)
        wasFactor<-TRUE
        x<-as.character(x)
      }
      else wasFactor<-FALSE
      if(is.character(x)){
        x[which(is.na(x))]<-as.character(newValue)
        if(wasFactor) x<-factor(x,levels=c(waslevels,as.character(newValue))) #keeps order of previous factors
      }
      else x[which(is.na(x))]<- newValue #assume numeric if not character/factor
      return(x)
    }
    #have to do this; otherwise makes them all characters if use apply...
    cnames<-colnames(sData)
    sData<-do.call("data.frame",lapply(seq_len(ncol(sData)),function(ii){fixNAFunction(sData[,ii],newValue=newValue)}))
    colnames(sData)<-cnames
    
  }
  return(sData)
}

.unnameClusterSlots<-function(ce){
  names(ce@clusterLegend)<-names(ce@clusterInfo)<-names(ce@clusterTypes)<-NULL
  return(ce)
}

#.convertToNum works on vector and is basically as.numeric, so that if character valued numbers, will get them, but other wise, arbitrary values based on making input x a factor (so alphabetical).
.convertToNum<-function(x){
  nms<-names(x)
  if(is.factor(x)){
    x<-as.character(x)
  }
  # FIXME what is this doing???
  if(is.character(x)){
    op <- options(warn=2)
    #if character values convert to numeric, will use that, otherwise convert to factor first
    test <- try(as.numeric(x) ,silent=TRUE)
    if(inherits(test,"try-error")) x<-as.numeric(factor(x))
    else x<-test
    options(op)
    
  }
  names(x)<-nms
  return(x)
}
##Universal way to convert matrix of clusters (of any value) into integers, preserving -1, -2 values. Makes them adjacent values. NAs are converted into -1 values.
.makeIntegerClusters<-function(clMat){
  if(!is.matrix(clMat)) stop("must give matrix input") #means must be either character or numeric values.
  fun<-function(x){ #make it numbers from 1:length(x), except for -1,-2
    #id special values of -1,-2
    isChar<-as.character(is.character(x))
		
		#make NA values NA
		if(any(is.na(x))){
			if(is.character(x)) x[is.na(x)]<-"-1" 
			else x[is.na(x)]<- -1 
		}
    wh1<-switch(isChar,"TRUE"= x =="-1","FALSE"= x ==-1)
    wh2<-switch(isChar,"TRUE"= x =="-2","FALSE"= x ==-2)
    wh<-wh1 | wh2
    vals<-sort(unique(x[!wh])) #if input was 1:n will hopefully make output 1:n matching
    y<-match(x,vals)
    y[wh1]<- -1
    y[wh2]<- -2
    return(y)
  }
  
  if(!is.null(dim(clMat)) && ncol(clMat)>1){
    x<-apply(clMat,2,fun)
    if(is.null(dim(x))) x<-matrix(x,nrow=1) #in case clMat was matrix with 1 row
    return(x  )
  } 
  else{
    if(is.matrix(clMat)) clMat<-clMat[,1]
    return(matrix(fun(clMat),ncol=1))
  }
}
## Universal way to convert matrix of clusters into default colorLegend
## returns  list(colorList=colorList,numClusters=clMat,facClusters=clMat))
## If giving clNumMat, then will not check for consecutive integers, but will use entries in clNumMat as the clusterIds
## If matchClusterLegend given will use names and colors of matchClusterLegend
## if matchTo="numIds", assume that 'clusterIds' column matches clNumMat (so only makes sense if giving clNumMat); will keep clusterId and color of clusterLegend -- NOT name
## if matchTo="givenIds", assume 'clusterIds' matches entries of clMat; will keep name and color of clusterLegend
## if matchTo="name" assume "name" matches entries of clMat; will keep name and color of clusterLegend
## (will only match if finds name or clusterIds in matchClusterLegend, depending on value of matchTo)
## Assumes that clMat has already had continuous columns removed
## Note, does drop levels, so returned datasets as well as color legend doesn't include missing factor levels
.makeColors<-function(clMat, colors,clNumMat=NULL,unassignedColor="white",missingColor="grey", distinctColors=FALSE,
                      matchClusterLegend=NULL,matchTo=c("givenIds","numIds","name")){ 
  
  matchTo<-match.arg(matchTo)
  if(!is.null(matchClusterLegend)){
    if(!is.list(matchClusterLegend) ) matchClusterLegend<-.convertToClusterLegend(matchClusterLegend)
    if("matchTo"%in% c("numIds","givenIds")) reqNames<-c("color","clusterIds") else reqNames<-c("color","name")
    ch<-.checkClusterLegendList(matchClusterLegend,allowNames=TRUE,reqNames=reqNames)
    if(!is.logical(ch)){
      #try again in case in aheatmap format
      if(!all(sapply(matchClusterLegend, function(x) {!is.null(dim(x))}))) matchClusterLegend<-.convertToClusterLegend(matchClusterLegend)
      ch<-.checkClusterLegendList(matchClusterLegend,allowNames=TRUE,reqNames=reqNames)
      if(!is.logical(ch)) stop(ch)
      
    }
		if(is.null(names(matchClusterLegend)) && !is.null(colnames(clMat)) && length(matchClusterLegend)==ncol(clMat)){
			names(matchClusterLegend)<-colnames(clMat)
		}
  }
  if(ncol(clMat)==1) distinctColors<-FALSE
  if(any(apply(clMat,2,function(x){length(unique(x))})>length(colors))) warning("too many clusters to have unique color assignments")
    ###(not sure why this simpler code doesn't give back data.frame with factors: annCol<-apply(clMat,2,function(x){factor(x)}))
		
	cNames<-colnames(clMat)
	#if give clNumMat, can skip this step, save time
	if(is.null(clNumMat)) clNumMat<-.makeIntegerClusters(as.matrix(clMat)) #require to not skip a integer value
	allFactors<-all(apply(clMat,2,function(x){is.factor(x)}))
	if(!allFactors){
		clMat<-do.call("data.frame", lapply(seq_len(ncol(clMat)), function(ii){ 
		  x<-clMat[,ii]
		  if(!is.factor(x)) return(factor(clMat[,ii]))
		  else return(droplevels(x))
		  }))
  	names(clMat)<-cNames
	}
	
  maxPerCol<-apply(clNumMat,2,max) #max cluster value (not including -1,-2)
  currcolors<-rep(colors,length= max(c(length(colors),sum(maxPerCol)))) #make sure don't run out of colors. Put max in case none are >0
  if(distinctColors) maxPreviousColor<-c(0,head(cumsum(maxPerCol),-1)) 
		
	#make a clusterLegend list
	perColumnFunction<-function(ii){
    facInt<-clNumMat[,ii] #assumes adjacent numbers
    nVals<-max(c(facInt,0))
    facOrig<-clMat[,ii] #assumes factors
    matchName<-cNames[[ii]]
		if(!is.null(matchClusterLegend[[matchName]])){
	    colMat<-matchClusterLegend[[matchName]]
			if(matchTo=="numIds"){
			  m<-match(colMat[,"clusterIds"],as.character(facInt))
			  cols<-cbind(colMat[,c("clusterIds","color"),drop=FALSE],"name"=as.character(facOrig)[m])
			}
			else if(matchTo=="givenIds"){
			  m<-match(colMat[,"clusterIds"],as.character(facOrig))
			 cols<-cbind("clusterIds"=facInt[m],colMat[,c("color","name"),drop=FALSE]) 			}
			else{
			  m<-match(colMat[,"name"],as.character(facOrig))
			  cols<-cbind("clusterIds"=facInt[m],colMat[,c("color","name"),drop=FALSE])
			}
	    if(any(is.na(m))) cols<-cols[!is.na(m),,drop=FALSE] #incase given legend has names/Ids not found in data
	    #in case given legend doesn't have values found in data
	    whMissing<-which(!as.character(facInt)[facInt>0] %in% cols[,"clusterIds"])
	    if(length(whMissing)>0){ 
	      #remove colors already in colMat
	      colors<-currcolors
	      mColors<-.matchColors(colors,cols[,"color"])
	      if(any(!is.na(mColors))) colors<-colors[-!is.na(mColors)]
	      mAdd<-match(unique(facInt[facInt>0][whMissing]),facInt) #make unique value
	      addMat<-cbind("clusterIds"=as.character(facInt[mAdd]),"color"=colors[seq_along(mAdd)],"name"=as.character(facOrig[mAdd]))
        cols<-rbind(cols,addMat)
	     }
	    
		}
		else{
		  if(nVals>0){
		    if(distinctColors){
		      add<-maxPreviousColor[[ii]]
		      colors<-currcolors[seq_len( nVals)+add]
		    }
		    else colors<-currcolors[seq_len( nVals)]
        #make matrix for non -1 values
		    clIds<-sort(unique(facInt[facInt>0]))
        m<-match(clIds,facInt)
		    cols<-cbind(
		      "clusterIds"=as.character(clIds),
		      "color"=colors,
		      "name"=as.character(facOrig[m])
		    )
		  }
		  else{ #all <0 cluster values
		    cols<-matrix(nrow=0,ncol=3)
		    colnames(cols)<-c("clusterIds","color","name")
		  }
	    if(any(facInt== -1)) 
				cols<-rbind(cols,c("clusterIds"="-1","color"=unassignedColor,"name"="-1") )
	    if(any(facInt== -2)) 
				cols<-rbind(cols,c("clusterIds"="-2","color"=missingColor,"name"="-2") )
		}
    cols<-cols[order(cols[,"name"]),,drop=FALSE]
		rownames(cols)<-NULL
    return(cols)
	}
	colorList<-lapply(seq_len(ncol(clNumMat)),FUN=perColumnFunction)
  names(colorList)<-cNames
  return(list(colorList=colorList,numClusters=clNumMat,facClusters=clMat))
}

#' @importFrom grDevices col2rgb 
.matchColors<-function(x,y){
  #match x colors to y colors, only with rgb definitions
  if(!is.character(x) || !is.character(y)) stop("colors must be character values")
  xrgb<-col2rgb(x)
  yrgb<-col2rgb(y)
  match(data.frame(xrgb), data.frame(yrgb))
}






