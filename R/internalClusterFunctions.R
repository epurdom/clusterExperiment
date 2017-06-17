#######
#Internal algorithms for clustering
#######
#check what type
.checkAlgType<-function(clusterFunction){
	##These return lists of indices of clusters satisifying alpha criteria
	if(clusterFunction=="tight") type<-"01"
	if(clusterFunction=="hierarchical01") type<-"01"
	if(clusterFunction=="hierarchicalK") type<-"K"
	if(clusterFunction=="pam") type<-"K"
	return(type)
}
#convert list output into cluster vector.
.convertClusterListToVector<-function(clusterList,N)
{
    clust.id <- rep(-1, N)
    nfound<-length(clusterList)
    if(nfound>0){
        #make cluster ids in order of when found
        for (i in 1:length(clusterList)) clust.id[clusterList[[i]]] <- i 
    }
    return(clust.id)
}
#Note, only returns 'both' if inputType is not given...otherwise picks
.checkXDissInput<-function(x,diss,inputType=NA,algType,checkDiss=TRUE){
  if(is.null(x) & is.null(diss)) stop("must give either x or diss argument")
  #  if(!is.null(x) & !is.null(diss)) stop("cannot give both x and diss argument")
  if(!is.null(x) & is.null(diss)) input<-"X"
  if(!is.null(x) & !is.null(diss)) input<-"both"
  if(is.null(x) & !is.null(diss)) input<-"diss"
  if(input %in% c("diss","both") & checkDiss) .checkDissFunction(diss,algType=algType)
  if(input == "both" && ncol(x)!=ncol(diss)) stop("ncol(x)!=ncol(diss): if both x and diss given then must have compatible dimensions.") 
  if(!is.na(inputType)){
	  if(input=="both"){
		  if(inputType=="diss") input<-"diss"
		  if(inputType=="X") input<-"X"
		  if(inputType=="either") input<-"diss"	#if both given and both acceptable, use diss.	  
	  }
	  if(input == "diss" & inputType=="X") stop("given clusterFunction/classifyFuntion only takes a X matrix")
	  #commented this out, because actually want the ability to use distFunction to make default diss if missing one. 
	#  if(input == "X" & inputType=="diss") stop("given clusterFunction/classifyFuntion only takes dissimilarity matrix")
  	
  }
 
  return(input)
}
.makeDiss<-function(x,distFunction,algType,checkDiss){
  if(!is.function(distFunction)){
	  if(length(distFunction)>1) stop("if distFunction is not a function, it must be of length 1")
	  if(is.character(distFunction)){
		  distFunction<-get(distFunction,envir=globalenv())
	  }else if(is.na(distFunction)){
	      distFunction<-switch(algType, "01"=function(x){(1-cor(t(x)))/2}, "K"=function(x){dist(x)})
	  }else stop("if distFunction is not a function, it must be either NA or a character")
  } 
  D<-try(as.matrix(distFunction(t(x))))	#distances assumed to be of observations on rows
  if(inherits(D,"try-error")) stop("input distance function gives error when applied to x")
  if(!all(dim(D) == c(ncol(x),ncol(x)))) stop("input distance function must result in a ",ncol(x),"by",ncol(x),"matrix of distances")
  if(checkDiss) .checkDissFunction(D,algType=algType)
  return(D)
	  
	
}
.checkDissFunction<-function(D,algType){
	if(any(is.na(as.vector(D)))) stop("NA values found in dissimilarity matrix (could be from too small of subsampling if classifyMethod!='All', see documentation of subsampleClustering)")
	if(any(is.na(D) | is.nan(D) | is.infinite(D))) stop("Dissimilarity matrix contains either NAs, NANs or Infinite values.")
	if(any(D<0)) stop("Dissimilarity matrix must have strictly positive values")
	if(any(diag(D)!=0)) stop("Dissimilarity matrix must have zero values on the diagonal")
	if(!all(D==t(D))) stop("Dissimilarity matrix must result in a symmetric matrix")
	if(algType=="01" & any(D>1)) stop("distance function must give values between 0 and 1 which algorithm type of the ClusterFunction object is '01'")
}




.clusterVectorToList<-function(vec){
    clList<-tapply(1:length(vec),vec,function(x){x},simplify=FALSE)
	whNotAssign<-which(sapply(clList,function(x){all(vec[x]== -1)}))
	if(length(whNotAssign)>1) stop("Internal coding error in removing unclustered samples")
    if(length(whNotAssign)>0) clList<-clList[-whNotAssign]	
}
.clusterListToVector<-function(ll,N){
	if(length(ll)==0) return(rep(-1,N))
	else{
		names(ll)<-as.character(1:length(ll))
		clId<-lapply(1:length(ll),function(ii){rep(ii,length=length(ll[[ii]]))})
		clVal<-unlist(clId)
		clInd<-unlist(ll)
		clusterVec<-rep(-1,length=N)
		clusterVec[clInd]<-clVal
		return(clusterVec)
		
	}
	
}
.orderByAlpha<-function(res,S)
{
	if(length(res)>0){
		alphaMax<-unlist(lapply(res, function(x){
			vals<-lower.tri(S[x,x]) #don't grab diag
			1-min(vals) #max(alpha)=1-min(S)
		}))
	    res <- res[order(alphaMax, decreasing=TRUE)]

	}
	else return(res)
}

.makeDataArgs<-function(dataInput,funInput,xData,dissData){
	if(dataInput=="X"){
		if(funInput=="diss") stop("Internal coding error: should have caught that wrong data input ('X') for this clusterFunction")
  		argsClusterList<-switch(funInput,"X"=list(x=xData), "either"=list(diss=NULL,x=xData))	
	}
	if(dataInput=="diss"){
		if(funInput=="X") stop("Internal coding error: should have caught that wrong data input ('diss') for this clusterFunction")
  		argsClusterList<-switch(funInput,"diss"=list(diss=dissData), "either"=list(diss=dissData,x=NULL)	)	
	}
	return(argsClusterList)	
}