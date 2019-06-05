#######
#Internal algorithms for cluster functions
#######
#convert list output into cluster vector.
.convertClusterListToVector <- function(clusterList, N)
{
    clust.id <- rep(-1, N)
    nfound <- length(clusterList)
    if (nfound > 0) {
        #make cluster ids in order of when found
        for (i in seq_along(clusterList))
            clust.id[clusterList[[i]]] <- i
    }
    return(clust.id)
}


.checkDissFunction <- function(D, algType = NA) {
    if (anyNA(D))
        stop(
            "NA values found in dissimilarity matrix (could be from too small of subsampling if classifyMethod!='All', see documentation of subsampleClustering)"
        )
    if (any(is.nan(D) |
            is.infinite(D)))
        stop("Dissimilarity matrix contains either NAs, NANs or Infinite values.")
    if (any(D < 0))
        stop("Dissimilarity matrix must have strictly positive values")
    if (any(diag(D) != 0))
        stop("Dissimilarity matrix must have zero values on the diagonal")
    if (!all(D == t(D)))
        stop("Dissimilarity matrix must result in a symmetric matrix")
    if (algType == "01" &
        any(D > 1))
        stop(
            "distance function must give values between 0 and 1 when algorithm type of the ClusterFunction object is '01'"
        )
}




.makeDiss <- function(x, distFunction, algType, checkDiss) {
    if (!is.function(distFunction)) {
        if (length(distFunction) > 1)
            stop("if distFunction is not a function, it must be of length 1")
        if (is.character(distFunction)) {
            distFunction <- get(distFunction, envir = globalenv())
        } else if (is.na(distFunction)) {
            distFunction <-
                switch(
                    algType,
                    "01" = function(x) {
                        (1 - cor(t(x))) / 2
                    },
                    "K" = function(x) {
                        dist(x)
                    }
                )
        } else
            stop("if distFunction is not a function, it must be either NA or a character")
    }
    ###Add data.matrix here for HDF5, not optimized.
    D <-
        try(as.matrix(distFunction(data.matrix(t(x)))))
    #distances assumed to be of observations on rows
    if (inherits(D, "try-error"))
        stop("input distance function gives error when applied to x")
    if (!all(dim(D) == c(ncol(x), ncol(x))))
        stop(
            "input distance function must result in a ",
            ncol(x),
            "by",
            ncol(x),
            "matrix of distances"
        )
    if (checkDiss)
        .checkDissFunction(D, algType = algType)
    return(D)
    
    
}


.clusterVectorToList <- function(vec) {
    clList <- tapply(seq_along(vec), vec, function(x) {
        x
    }, simplify = FALSE)
    whNotAssign <- which(sapply(clList, function(x) {
        all(vec[x] == -1)
    }))
    if (length(whNotAssign) > 1)
        stop("Internal coding error in removing unclustered samples")
    if (length(whNotAssign) > 0)
        clList <- clList[-whNotAssign]
}
.clusterListToVector <- function(ll, N) {
    if (length(ll) == 0)
        return(rep(-1, N))
    else{
        names(ll) <- as.character(seq_along(ll))
        clId <-
            lapply(seq_along(ll), function(ii) {
                rep(ii, length = length(ll[[ii]]))
            })
        clVal <- unlist(clId)
        clInd <- unlist(ll)
        clusterVec <- rep(-1, length = N)
        clusterVec[clInd] <- clVal
        return(clusterVec)
        
    }
    
}
.orderByAlpha <- function(res, S)
{
    if (length(res) > 0) {
        alphaMax <- unlist(lapply(res, function(x) {
            vals <- lower.tri(S[x, x]) #don't grab diag
            1 - min(vals) #max(alpha)=1-min(S)
        }))
        res <- res[order(alphaMax, decreasing = TRUE)]
        
    }
    else
        return(res)
}


# if inputType is vector, returns all those that match given input
# if inputType is NA, just returns those that are found (i.e. doesn't narrow them down.)
.checkCFInput <- function(x,
                          diss,
                          cat,
                          inputType = NA,
                          algType,
                          checkDiss) {
    if (!all(inputType %in% .inputTypes))
        stop(sprintf(
            "not all inputType values match allowable values (%s)",
            paste(.inputTypes, ", ")
        ))
    if (is.null(x) &
        is.null(diss) &
        is.null(cat))
        stop("must give either x or diss or cat argument")
    #  if(!is.null(x) & !is.null(diss)) stop("cannot give both x and diss argument")
    input <- c()
    if (!is.null(x))
        input <- c(input, "X")
    if (!is.null(diss))
        input <- c(input, "diss")
    if (!is.null(cat))
        input <- c(input, "cat")
    
    ##If multiple inputs, check that they are compatible dimensions
    bothXDiss <- all(c("X", "diss") %in% input)
    if (bothXDiss %in% input)
        &&
        ncol(x) != ncol(diss)) stop("ncol(x)!=ncol(diss): if both x and diss given then must have compatible dimensions.")
bothCatDiss<-all(c("X","diss") %in% input)
if(bothCatDiss %in% input) && ncol(cat)!=ncol(diss)) stop("ncol(cat)!=ncol(diss): if both cat and diss given then must have compatible dimensions.")
bothCatX<-all(c("X","cat") %in% input)
if(bothCatX %in% input) && ncol(cat)!=ncol(X)) stop("ncol(cat)!=ncol(X): if both cat and diss given then must have compatible dimensions.")

if(!is.na(inputType)){
    intersect<-intersect(inputType,input)
    if(length(intersect)==0){
        stop(sprintf("given clusterFunction/classifyFunction does not take as input matrices of type %s",paste(input,collapse=" or "))
    }
}
else{
    intersect<-input
}
if(any(intersect %in% c("diss")) & checkDiss) .checkDissFunction(diss,algType=algType)
return(intersect)
}

#' @param dataInput one of "X", "diss" to indicate type of data
#' @param funInput one of "X","diss" to indicate type function expects
#' @param xData the "X" data
#' @param dissData the "diss" data
#' @return returns list of arguments of the data with the corrected input that can be combined in a do.call to the function
.makeDataArgs <- function(dataInput,
                          funInput,
                          xData,
                          dissData,
                          catData) {
    if (dataInput == "X") {
        if (funInput == "diss")
            stop(
                "Internal coding error: should have caught that wrong data input ('X') for this clusterFunction"
            )
        argsClusterList <-
            switch(
                funInput,
                "X" = list(x = xData),
                "either" = list(diss = NULL, x = xData)
            )
    }
    if (dataInput == "diss") {
        if (funInput %in% c("X", "cat"))
            stop(
                "Internal coding error: should have caught that wrong data input ('diss') for this clusterFunction"
            )
        argsClusterList <-
            switch(
                funInput,
                "diss" = list(diss = dissData),
                "either" = list(diss = dissData, x = NULL)
            )
    }
    if (dataInput == "cat") {
        if (funInput %in% c("X", "diss"))
            stop(
                "Internal coding error: should have caught that wrong data input ('cat') for this clusterFunction"
            )
        argsClusterList <-
            switch(
                funInput,
                "cat" = list(cat = catData),
                "either" = list(diss = dissData, x = NULL)
            )
    }
    return(argsClusterList)
}

###This function checks the mainClusterArgs and subsampleArgs to make sure makes sense with combination of sequential, subsample, x, and diss given by the user.
#' @return If there is error, returns a character string describing error, otherwise returns list with necessary information:
#' inputClusterD
#' mainClusterArgs
#' subsampleArgs
.checkSubsampleClusterDArgs <-
    function(x,
             diss,
             cat,
             subsample,
             sequential,
             mainClusterArgs,
             subsampleArgs,
             checkDiss,
						 allowMakeDiss=FALSE, #if FALSE, will give error if not match. 
             warn = TRUE) {
        ########
        #checks for mainClustering stuff
        ########
				makeDiss<-FALSE
				if (!is.null(mainClusterArgs)) {
            if ("clusterFunction" %in% names(mainClusterArgs)) {
                #get clusterFunction for cluster D
                clusterFunction <- mainClusterArgs[["clusterFunction"]]
                if (is.character(clusterFunction))
                    clusterFunction <- getBuiltInFunction(clusterFunction)
                
                #Following input commands will return only X or Diss because gave the inputType argument...
                input <-
                    .checkCFInput(
                        x,
                        diss,
                        inputType = inputType(clusterFunction),
                        algType = algorithmType(clusterFunction),
                        checkDiss = checkDiss
                    )
                algType <- algorithmType(clusterFunction)
            }
            else{
                return(
                    "Must provide 'clusterFunction' for the mainClustering step to be able to run (give 'clusterFunction' argument via 'mainClusterArgs')"
                )
            }
            #-----
            #this check is done in mainClustering: removes certain options based on others
            # but want to do this check before run subsampling...also can give message that clearer that it refers to mainClustering.
            #-----
            reqArgs <- requiredArgs(clusterFunction)
            # remove 'k' is sequential is given
            if (sequential &
                length(reqArgs) > 0)
                reqArgs <- reqArgs[-which(reqArgs == "k")]
            # remove 'k' if choose 'findBestK=TRUE'
            if (length(reqArgs) > 0 &
                algorithmType(clusterFunction) == "K" &
                "findBestK" %in% names(mainClusterArgs)) {
                if (mainClusterArgs[["findBestK"]])
                    reqArgs <- reqArgs[-which(reqArgs == "k")]
            }
            
            #Check minSize valid
            if (!is.numeric(mainClusterArgs[["minSize"]]) ||
                mainClusterArgs[["minSize"]] < 0)
                return(
                    "Invalid value for the 'minSize' parameter in determining the minimum number of samples required in a cluster (in main clustering step)."
                )
            else
                mainClusterArgs[["minSize"]] <-
                    round(mainClusterArgs[["minSize"]]) #incase not integer.
            
            ## Check post-processing arguments
            doKPostProcess <- FALSE
            whPostProcessArgs <-
                which(names(mainClusterArgs) %in% getPostProcessingArgs(clusterFunction))
            if (length(whPostProcessArgs) > 0) {
							  ##FIXME: Do I need both checkArgs and warn?
                #get rid of wrong args passed because of user confusion between the two
                if (!is.null(mainClusterArgs[["postProcessNames"]]) && length(whPostProcessArgs) != length(mainClusterArgs[["postProcessNames"]]) & checkArgs & warn)
                    warning(
                        "Some arguments passed via '...' in mainClusterArgs in mainClustering step do not match the algorithmType of the given ClusterFunction object"
                    )
                postProcessArgs <- mainClusterArgs[whRightArgs]
                if (clusterFunction@algorithmType == "K") {
                    if ("findBestK"  %in% names(postProcessArgs) &&
                        postProcessArgs[["findBestK"]])
                        doKPostProcess <- TRUE
                    if ("removeSil"  %in% names(postProcessArgs) &&
                        postProcessArgs[["removeSil"]])
                        doKPostProcess <- TRUE
                }
                mainClusterArgs[["postProcessArgs"]] <- postProcessArgs
            }
						else mainClusterArgs[["postProcessArgs"]]<-NULL
            mainClusterArgs[["doKPostProcess"]] <- doKPostProcess
            
            #------
            # Check have required args for mainClustering
            #------
            if (length(reqArgs) > 0) {
                if (("clusterArgs" %in% names(mainClusterArgs) &
                     !all(reqArgs %in% names(mainClusterArgs[["clusterArgs"]]))) ||
                    !("clusterArgs" %in% names(mainClusterArgs)))
                    return(
                        paste(
                            "For the clusterFunction algorithm type ('",
                            algorithmType(clusterFunction),
                            "') given in 'mainClusterArgs', must supply arguments:",
                            reqArgs,
                            "These must be supplied as elements of the list of 'clusterArgs' given in 'mainClusterArgs'"
                        )
                    )
            }
            
            #------
            # Check if should calculate diss
            #------
            if (input == "X" &
                (clusterFunction@inputType == "diss" || doKPostProcess)) {
                if(allowMakeDiss){
									makeDiss <- TRUE
									if(warn) warning(sprintf("The following arguments to the mainClustering step require calculation of the n x n dissimilarity matrix: input=%s, clusterFunction@inputType=%s", input, clusterFunction@inputType))
									
								}
								else return("The clusterFunction requested in the mainClustering step requires a dissimilarity matrix as input")
                if (clusterFunction@inputType == "diss")
                    input <- "diss"
            }
        }
        
        
        ########
        # Adapt MainClusteringArgs if using Sequential
        ########
        if (sequential & !is.null(mainClusterArgs)) {
            # Remove argument 'k'
            # Reason: if subsample=FALSE, then need to change k of the mainClustering step for sequential. If subsample=TRUE, similarly set the k of mainClustering step to match that used in subsample. Either way, can't be predefined by user
            if ("clusterArgs" %in% names(mainClusterArgs)) {
                if ("k" %in% names(mainClusterArgs[["clusterArgs"]])) {
                    #remove predefined versions of k from both.
                    whK <- which(names(mainClusterArgs[["clusterArgs"]]) == "k")
                    if (warn)
                        warning(
                            "Setting 'k' in mainClusterArgs when sequential clustering is requested will have no effect."
                        )
                    mainClusterArgs[["clusterArgs"]] <-
                        mainClusterArgs[["clusterArgs"]][-whK]
                }
            } else{
                mainClusterArgs[["clusterArgs"]] <-
                    NULL #make it exist... not sure if I need this.
            }
        }
        #########
        # Checks related to subsample=TRUE
        #########
        if (subsample) {
            # Check that mainClustering cluster function takes input that is Diss
            # Reason: if subsampling, then the D from subsampling sent to the clusterFunction.
            if (!is.null(mainClusterArgs) && inputType(clusterFunction) == "X" )
                return(
                    "If choosing subsample=TRUE, the clusterFunction used in the mainClustering step must take input that is dissimilarity."
                )
						
            #--------
						# Evaluate subsample clusterFunction
						#--------
            if ("clusterFunction" %in% names(subsampleArgs)) {
                #---
                #Checks for when cluster function set by user for subsampling..
                #---
                subsampleCF <- subsampleArgs[["clusterFunction"]]
                if (is.character(subsampleCF))
                    subsampleCF <- getBuiltInFunction(subsampleCF)
                subsampleAlgType <- algorithmType(subsampleCF)
                
						    inputSub<-.checkCFInput(x=x, diss=diss, cat=cat,inputType=subsampleCF@inputType, checkDiss=checkDiss)
								
                # Check that subsample clustering function is of type 'k' if sequential=TRUE
                # Reason: seqCluster requires subsampling cluster function to be of type "K"
                if (sequential & algorithmType(subsampleCF) != "K") {
                    if (warn)
                        warning(
                            sprintf(
                                "If subsample=TRUE, sequentical clustering can only be implemented with a clusterFunction for subsampling that has algorithmType 'K'. See documentation of seqCluster. Will ignore this argument of subsampleArgs and set to default of %s",
                                default
                            )
                        )
                    subsampleArgs[["clusterFunction"]] <- default
                    diffSubsampleCF <- TRUE
                }
            }
            else if(!is.null(mainClusterArgs)){
	            #default clusterFunction for subsampled data, if not given.
	            default <- switch(inputSub,
	                              "X" = "kmeans",
	                              cat = "hier01",
	                              diss = "pam")
							
                #---
                #Checks for when no cluster function set for subsampling..
                #---
                if (!sequential || algorithmType(clusterFunction) == "K") {
                    mess <-
                        "a clusterFunction was not set for subsampleClustering -- set to be the same as the mainClustering step"
                    if (is.character(mainClusterArgs[["clusterFunction"]]))
                        mess <- sprintf("%s (%s)", mess, mainClusterArgs[["clusterFunction"]])
                    if (warn)
                        .mynote(mess)
                    subsampleArgs[["clusterFunction"]] <- clusterFunction
                    inputSubsample <- inputSub
                    diffSubsampleCF <- FALSE
                }
                else{
                    if (warn)
                        .mynote(
                            sprintf(
                                "a clusterFunction was not set for subsampleClustering (and sequential=TRUE means that it must be of type 'K' so cannot be set to that of mainClustering step). The clusterFunction for subsampling was set to the default of %s",
                                default
                            )
                        )
                    subsampleArgs[["clusterFunction"]] <- default
                    diffSubsampleCF <- TRUE
                }
                #Update for any changes from above
                subsampleCF <- subsampleArgs[["clusterFunction"]]
                if (is.character(subsampleCF))
                    subsampleCF <- getBuiltInFunction(subsampleCF)
                # Note: this will check diss again if checkDiss=TRUE
                # Reason: if algorithm on main is 01 and other isn't, need to check diss again.
                if (diffSubsampleCF) {
                    inputSubsample <-
                        .checkCFInput(
                            x,
                            diss,
                            inputType = inputType(subsampleCF),
                            algType = algorithmType(subsampleCF),
                            checkDiss = checkDiss
                        )
                }
                
            }
						else return("Must give 'clusterFunction' value to 'subsampleClustering'.")
            
							#Check the classify function
						if (is.null(subsampleCF@classifyFUN)) {
                if ("classifyMethod" %in% names(subsampleArgs) &&
                    subsampleArgs[["classifyMethod"]] != "InSample")
                    stop(
                        "Cannot set 'classifyMethod' to anything but 'InSample' if do not specify a clusterFunction in subsampleArgs that has a non-null classifyFUN slot"
                    )
                subsampleArgs[["classifyMethod"]] <- "InSample"
            }
						inputClassify<-.checkXDissInput(x, diss, cat inputType=subsampleCF@inputClassifyType, checkDiss=FALSE) #don't need to check it twice, even if asked for it!
						if(subsampleArgs[["classifyMethod"]]!="InSample"){
							if(inputClassify=="X" && subsampleCF@inputClassifyType=="diss")
								if(warn) warning(sprintf("The following arguments regarding the clusterFunction for subsampleArgs require calculation of the n x n dissimilarity matrix: classifyMethod=%s, inputClassify=%s, inputClassifyType=%s. Changing classifyMethod to 'InSample'", subsampleArgs[["classifyMethod"]], inputClassify, subsampleCF@inputClassifyType)
					    }
							subsampleArgs[["classifyMethod"]]<-"InSample"
				    }
					    
            #Reason: check subsampleArgs has required arguments for function, repeated from subsamplingClustering, but want it here before do calculations... if not, see if can borrow from mainClusterArgs
            ##------
            ##Check have required args for subsample. If missing, 'borrow' those args from mainClusterArgs.
            ##------
            reqSubArgs <- requiredArgs(subsampleCF)
            
            #Reason: sequential sets k for the subsampling via k0
            if (sequential & length(reqSubArgs) > 0)
                reqSubArgs <- reqSubArgs[-which(reqSubArgs == "k")]
            if (length(reqSubArgs) > 0)}
						 		if(!is.null(mainClusterArgs)) {
                #check if can borrow...
                if ("clusterArgs" %in% names(mainClusterArgs)) {
                    mainReqArgs <- requiredArgs(clusterFunction)
                    mainReqArgs <-
                        mainReqArgs[mainReqArgs %in% names(mainClusterArgs[["clusterArgs"]])]
                    if (!is.null(subsampleArgs) &&
                        "clusterArgs" %in% names(subsampleArgs)) {
                        #check if existing clusterArgs has required names already
                        #if not, give them those of mainClustering if exist.
                        if (!all(reqSubArgs %in% names(subsampleArgs[["clusterArgs"]]))) {
                            missingArgs <-
                                reqSubArgs[!reqSubArgs %in% names(subsampleArgs[["clusterArgs"]])]
                            missingArgs <- missingArgs[missingArgs %in% mainReqArgs]
                            
                        }
                        else
                            missingArgs <- c()
                    }
                    else{
                        missingArgs <- reqSubArgs[reqSubArgs %in% mainReqArgs]
                    }
                    if (length(missingArgs) > 0) {
                        subsampleArgs[["clusterArgs"]][missingArgs] <-
                            mainClusterArgs[["clusterArgs"]][missingArgs]
                        if (warn)
                            warning(
                                "missing arguments ",
                                missingArgs,
                                " provided from those in 'mainClusterArgs'"
                            )
                    }
                }
                #now check if got everything needed...
                if (("clusterArgs" %in% names(subsampleArgs) &
                     !all(reqSubArgs %in% names(subsampleArgs[["clusterArgs"]]))) ||
                    !("clusterArgs" %in% names(subsampleArgs)))
                    return(
                        paste(
                            "For the clusterFunction algorithm type ('",
                            algorithmType(subsampleCF),
                            "') given in 'subsampleArgs', must supply arguments:",
                            reqSubArgs,
                            ". These must be supplied as elements of the list of 'clusterArgs' given in 'subsampleArgs'"
                        )
                    )
									}
									else{
								    #-----
								    # Other Checks
								    #-----
								    
								    if(!all(reqSubArgs %in% names(subsampleArgs[["clusterArgs"]]))) return(paste("For this clusterFunction algorithm type ('",algorithmType(subsampleCF),"') must supply arguments",reqArgs,"as elements of the list of 'clusterArgs' in the 'subsampleArgs' argument."))
            }
						
					    
						
						}
            #Reason, if subsample=TRUE, user can't set distance function because use diss from subsampling.
            if(!is.null(mainClusterArgs) &&
                "distFunction" %in% names(mainClusterArgs) &&
                !is.na(mainClusterArgs[["distFunction"]])) {
                if (warn)
                    warning(
                        "if 'subsample=TRUE', 'distFunction' argument in mainClusterArgs is ignored."
                    )
                mainClusterArgs[["distFunction"]] <- NA
            }
            if (sequential) {
                #Reason: if subsampling, sequential goes over different k values, so user can't set k
                if ("clusterArgs" %in% names(subsampleArgs) &&
                    "k" %in% names(subsampleArgs[["clusterArgs"]])) {
                    #remove predefined versions of k from both.
                    whK <- which(names(subsampleArgs[["clusterArgs"]]) == "k")
                    if (warn)
                        warning(
                            "Setting 'k' in subsampleArgs when sequential=TRUE is called will have no effect."
                        )
                    subsampleArgs[["clusterArgs"]] <-
                        subsampleArgs[["clusterArgs"]][-whK]
                }
                if (!"clusterArgs" %in% names(subsampleArgs)) {
                    subsampleArgs[["clusterArgs"]] <- list() #make it if doesn't exist
                }
            }
        }
        else{
            #not subsample, additional 
            if (sequential & !is.null(mainClusterArgs)) {
                #Reason: if subsample=FALSE, and sequential then need to adjust K in the mainClustering step, so need algorithm of type K
                if (algorithmType(clusterFunction) != "K") {
                    return(
                        "if subsample=FALSE, sequentical clustering can only be implemented with a clusterFunction with algorithmType 'K'. See documentation of seqCluster."
                    )
                    subsampleArgs <-
                        subsampleArgs[-which(names(subsampleArgs) == "clusterFunction")]
                }
                #Reason: subsample=FALSE can't do sequential clustering and findBestK=TRUE because need to remove cluster based on testing many k and finding stable, and if not doing it over subsample, then do it over actual clustering
                if (algorithmType(clusterFunction) == "K") {
                    if ("findBestK" %in% names(mainClusterArgs)) {
                        if (mainClusterArgs[["findBestK"]])
                            return(
                                "Cannot do sequential clustering where subsample=FALSE and 'findBestK=TRUE' is passed to the mainClustering step via mainClusterArgs. See help documentation of seqCluster."
                            )
                    }
                }
            }
            
        }
        
        return(
            list(
                inputClusterD = input,
                mainClusterArgs = mainClusterArgs,
                subsampleArgs = subsampleArgs,
								makeDiss=makeDiss
            )
        )
        
    }
