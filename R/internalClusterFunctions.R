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
.checkDissFunction<-function(D,algType=NA){
	if(any(is.na(as.vector(D)))) stop("NA values found in dissimilarity matrix (could be from too small of subsampling if classifyMethod!='All', see documentation of subsampleClustering)")
	if(any(is.na(D) | is.nan(D) | is.infinite(D))) stop("Dissimilarity matrix contains either NAs, NANs or Infinite values.")
	if(any(D<0)) stop("Dissimilarity matrix must have strictly positive values")
	if(any(diag(D)!=0)) stop("Dissimilarity matrix must have zero values on the diagonal")
	if(!all(D==t(D))) stop("Dissimilarity matrix must result in a symmetric matrix")
	if(algType=="01" & any(D>1)) stop("distance function must give values between 0 and 1 when algorithm type of the ClusterFunction object is '01'")
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

###This function checks the clusterDArgs and subsampleArgs to make sure make sense with combination of sequential, subsample, x, and diss given by the user. If error, returns a character string describing error, otherwise returns list with necessary information.
.checkSubsampleClusterDArgs<-function(x,diss,subsample,sequential,clusterDArgs,subsampleArgs,checkDiss){
    ########
 	#checks for clusterD stuff
    ########
     if("clusterFunction" %in% names(clusterDArgs)){
      #get clusterFunction for cluster D
  	  clusterFunction<-clusterDArgs[["clusterFunction"]]
  	  if(is.character(clusterFunction)) clusterFunction<-getBuiltInClusterFunction(clusterFunction)
  	  
	  #Following input commands will return only X or Diss because gave the inputType argument...
  	  input<-.checkXDissInput(x, diss, inputType=inputType(clusterFunction), algType=algorithmType(clusterFunction), checkDiss=checkDiss)
	  algType<-algorithmType(clusterFunction)
    }
    else{
  	  return("Must provide 'clusterFunction' for the clusterD step to be able to run (give 'clusterFunction' argument via 'clusterDArgs')")
    }
	#this check is done in clusterD, but want to do this check before run subsampling...also can give message that clearer that it refers to clusterD.
 	reqArgs<-requiredArgs(clusterFunction)
	if(sequential & length(reqArgs)>0) reqArgs<- reqArgs[-which(reqArgs=="k")]
	#remove required args not needed if certain postProcessArgs are given:
	# don't need define 'k' if choose 'findBestK=TRUE'
	if(length(reqArgs)>0 & algorithmType(clusterFunction)=="K" & "findBestK" %in% names(clusterDArgs)){
		if(clusterDArgs[["findBestK"]]) reqArgs<-reqArgs[-which(reqArgs=="k")]
	}
	if(length(reqArgs)>0){
		if(("clusterArgs"%in% names(clusterDArgs) & !all(reqArgs %in% names(clusterDArgs[["clusterArgs"]]))) || !("clusterArgs"%in% names(clusterDArgs))) return(paste("For the clusterFunction algorithm type ('",algorithmType(clusterFunction),"') given in 'clusterDArgs', must supply arguments:",reqArgs,"These must be supplied as elements of the list of 'clusterArgs' given in 'clusterDArgs'"))
	} 
    if(sequential){
		#Reason, if subsample=FALSE, then need to change k of the clusterD step for sequential. If subsample=TRUE, similarly set the k of clusterD step to match that used in subsample. Either way, can't be predefined by user
    	  if("clusterArgs" %in% names(clusterDArgs)){
			  if("k" %in% names(clusterDArgs[["clusterArgs"]]) ){
          	#remove predefined versions of k from both.
          	whK<-which(names(clusterDArgs[["clusterArgs"]])=="k")
          	warning("Setting 'k' in clusterDArgs when sequential clustering is requested will have no effect.")
          	clusterDArgs[["clusterArgs"]]<-clusterDArgs[["clusterArgs"]][-whK]
    			}
		}else{
			clusterDArgs[["clusterArgs"]]<-NULL #make it exist... not sure if I need this.
		}
	}
	#########
	# Checks related to subsample=TRUE
	#########
    if(subsample){
    	#Reason: if subsampling, then the D from subsampling sent to the clusterFunction.
    	if(inputType(clusterFunction)=="X") return("If choosing subsample=TRUE, the clusterFunction used in the clusterD step must take input that is dissimilarity.")
     	if("clusterFunction" %in% names(subsampleArgs)){	    
		  subsampleCF<-subsampleArgs[["clusterFunction"]]
		  if(is.character(subsampleCF)) subsampleCF<-getBuiltInClusterFunction(subsampleCF)
		  subsampleAlgType<-algorithmType(subsampleCF)
		  #Reason: seqCluster requires subsampling cluster function to be of type "K"
		  if(sequential & algorithmType(subsampleCF)!="K"){
			  warning("If subsample=TRUE, sequentical clustering can only be implemented with a clusterFunction for subsampling that has algorithmType 'K'. See documentation of seqCluster. Will ignore this argument of subsampleArgs and use the default clustering function.")
			  subsampleArgs<-subsampleArgs[-which(names(subsampleArgs)=="clusterFunction")]
		  }
  		  inputSubsample<-.checkXDissInput(x,diss, inputType=inputType(subsampleCF),  algType=algorithmType(subsampleCF), checkDiss=checkDiss) #if algorithm on one is 01 and other isn't, need to check diss again.
  		  diffSubsampleCF<-TRUE
		}
#  		else stop("must provide clusterFunction to subsampleArgs if subsample=TRUE")
		#this makes default to be same as clusterD
	  	else{
			warning("a clusterFunction was not set for subsampleClustering -- set to be the same as the clusterD step.")
			subsampleArgs[["clusterFunction"]]<-clusterFunction
	  		subsampleCF<-clusterFunction
	  		inputSubsample<-input
	  		diffSubsampleCF<-FALSE
	  	}
		if(is.null(subsampleCF@classifyFUN)){
			if("classifyMethod" %in% names(subsampleArgs) && subsampleArgs[["classifyMethod"]]!="InSample") stop("Cannot set 'classifyMethod' to anything but 'InSample' if do not specify a clusterFunction in subsampleArgs that has a non-null classifyFUN slot")
			subsampleArgs[["classifyMethod"]]<-"InSample"
		}
	  #Reason: check subsampleArgs has required arguments for function, repeated from subsamplingClustering, but want it here before do calculations... if not, see if can borrow from clusterDArgs
		##------
		##Check have required args for subsample. If missing, 'borrow' those args from clusterDArgs.
		##------
		reqSubArgs<-requiredArgs(subsampleCF)
  	
		#Reason: sequential sets k for the subsampling via k0
		if(sequential & length(reqSubArgs)>0) reqSubArgs<- reqSubArgs[-which(reqSubArgs=="k")]
		if(length(reqSubArgs)>0){
			#check if can borrow...
			if("clusterArgs" %in% names(clusterDArgs)){
				clusterDReqArgs<-requiredArgs(clusterFunction)
				clusterDReqArgs<-clusterDReqArgs[clusterDReqArgs%in%names(clusterDArgs[["clusterArgs"]])]
				if(!is.null(subsampleArgs) && "clusterArgs" %in% names(subsampleArgs)){
					#check if existing clusterArgs has required names already
					#if not, give them those of clusterD if exist.
					if(!all(reqSubArgs %in% names(subsampleArgs[["clusterArgs"]]))) {
						missingArgs<-reqSubArgs[!reqSubArgs%in%names(subsampleArgs[["clusterArgs"]])]
						missingArgs<-missingArgs[missingArgs%in%clusterDReqArgs]
						
				    }
					else missingArgs<-c()
				} 	
				else{
					missingArgs<-reqSubArgs[reqSubArgs%in%clusterDReqArgs]
				}
				if(length(missingArgs)>0){
					subsampleArgs[["clusterArgs"]][missingArgs]<-clusterDArgs[["clusterArgs"]][missingArgs]
					warning("missing arguments ",missingArgs," provided from those in 'clusterDArgs'")
				}
			}
			#now check if got everything needed...
			if(("clusterArgs" %in% names(subsampleArgs) & !all(reqSubArgs %in% names(subsampleArgs[["clusterArgs"]]))) || !("clusterArgs"%in% names(subsampleArgs))) 					return(paste("For the clusterFunction algorithm type ('",algorithmType(subsampleCF),"') given in 'subsampleArgs', must supply arguments:",reqSubArgs,". These must be supplied as elements of the list of 'clusterArgs' given in 'subsampleArgs'"))	
	  	  }
		  #Reason, if subsample=TRUE, user can't set distance function because use diss from subsampling.
		    if(!is.null(clusterDArgs) && "distFunction" %in% names(clusterDArgs) && !is.na(clusterDArgs[["distFunction"]])){
		        warning("if 'subsample=TRUE', 'distFunction' argument in clusterDArgs is ignored.")
		        clusterDArgs[["distFunction"]]<-NA
		    }
			if(sequential){
			  #Reason: if subsampling, sequential goes over different k values, so user can't set k
			   if("clusterArgs" %in% names(subsampleArgs) && "k" %in% names(subsampleArgs[["clusterArgs"]])){
		  	      #remove predefined versions of k from both.
		  	      whK<-which(names(subsampleArgs[["clusterArgs"]])=="k")
		  	      warning("Setting 'k' in subsampleArgs when sequential=TRUE is called will have no effect.")
		  	      subsampleArgs[["clusterArgs"]]<-subsampleArgs[["clusterArgs"]][-whK]
			  	}
		  	  if(!"clusterArgs" %in% names(subsampleArgs) ){
		  	  	  subsampleArgs[["clusterArgs"]]<-list() #make it if doesn't exist
		  	  }				 
		}
	}
	else{ #not subsample
	    if(sequential){
			#Reason: if subsample=FALSE, and sequential then need to adjust K in the clusterD step, so need algorithm of type K
	    	  if(algorithmType(clusterFunction) != "K"){
	    		  return("if subsample=FALSE, sequentical clustering can only be implemented with a clusterFunction with algorithmType 'K'. See documentation of seqCluster.")
	    		  subsampleArgs<-subsampleArgs[-which(names(subsampleArgs)=="clusterFunction")]	
	    	  }
		    #Reason: subsample=FALSE can't do sequential clustering and findBestK=TRUE because need to remove cluster based on testing many k and finding stable, and if not doing it over subsample, then do it over actual clustering
		    if(algorithmType(clusterFunction) == "K"){
		      if("findBestK" %in% names(clusterDArgs)){
		        if(clusterDArgs[["findBestK"]]) return("Cannot do sequential clustering where subsample=FALSE and 'findBestK=TRUE' is passed to the clusterD step via clusterDArgs. See help documentation of seqCluster.")
		      }
		    }
	    }

	}

    return(list(inputClusterD=input,clusterDArgs=clusterDArgs,subsampleArgs=subsampleArgs))
  
}


