#######
#Internal algorithms for cluster functions
#######

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
    return(clList)
}
.clusterListToVector <- function(clusterList, N)
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



.defaultDiss<-c("euclidean", "maximum", "manhattan", "canberra", "binary" ,"minkowski")
.myDefaultDiss<-c("01","K")
.diss01<-function(x) {
    (1 - cor(t(x))) / 2
}
.dissK<-function(x) {
    dist(x) #distance finds between rows of matrix
}

## x is a nfeatures x nsamples matrix
.makeDiss <- function(x, distFunction, algType, checkDiss) {
    if (!is.function(distFunction)) {
        if (length(distFunction) > 1)
            stop("if distFunction is not a function, it must be of length 1")
        if (is.character(distFunction) & distFunction!="default") {
            ## Handle defaults in `dist` function
            if(distFunction %in% .defaultDiss){
                distMethod<-distFunction
                distFunction<-function(x){dist(x,method=distMethod)}
            }
            else{
                if(distFunction %in% .myDefaultDiss){
                    if(distFunction=="01") distFunction<-.diss01
                    else if(distFunction=="K") distFunction<-.dissK
                }
                else distFunction <- get(distFunction, envir = globalenv())
            }
        } else if (is.na(distFunction) || distFunction=="default") {
            distFunction <-
                switch(
                    algType,
                    "01" = .diss01,
                    "K" = .dissK
                )
        } else
            stop("if distFunction is not a function, it must be either NA or a character")
    }
    ### FIXME: HDF5 optimization could be improved?
    ### Add data.matrix here for HDF5, not optimized.
    D <-try(as.matrix(distFunction(data.matrix(t(x)))))
    #distances assumed to be of observations on rows
    if (inherits(D, "try-error"))
        stop("given distance function gives error when applied to x")
    if (!all(dim(D) == c(ncol(x), ncol(x))))
        stop(
            "given distance function must result in a ",
            ncol(x),
            "by",
            ncol(x),
            "matrix of distances"
        )
    if (checkDiss)
        .checkDissFunction(D, algType = algType)
    return(D)
    
    
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


# if inputType is character value, returns if matches desiredInputType, otherwise throws appropriate errors. 

.checkCFInput <- function(inputType, desiredInputType) {
    if (!all(desiredInputType %in% .inputTypes))
        stop(sprintf(
            "not all values in 'desiredInputType' match allowable values (%s)",
            paste(.inputTypes, ", ")
        ))
    intersect<-intersect(desiredInputType,inputType)
    if(length(intersect)==0){
        stop(sprintf("given clusterFunction/classifyFunction does not take as input matrices of type %s",paste(inputType,collapse=" or ")))
    }
    return(inputType)
}

###This function checks the mainClusterArgs and subsampleArgs to make sure makes sense with combination of sequential, subsample, x, and diss given by the user.
#' @return If there is error, returns a character string describing error, otherwise returns list with necessary information:
#' mainClusterArgs
#' subsampleArgs
#' makeDiss (if allowMakeDiss=FALSE, will always be FALSE)
#' Note that warnings will be printed if warn=TRUE (different than errors...)
.checkArgs <-
    function(inputType,
             main=TRUE,
             subsample,
             sequential,
             mainClusterArgs,
             subsampleArgs,
             seqArgs=NULL,
             allowMakeDiss=FALSE,
             warn = TRUE) {
        ########
        #checks for mainClustering stuff
        ########
        inputType<-match.arg(inputType,.inputTypes)
        doDiss<-FALSE #If returns TRUE, then was a mismatch, but could be fixed by X -> Diss
        doDissPost<-FALSE
        if (main) {
            if ("clusterFunction" %in% names(mainClusterArgs)) {
                #get clusterFunction for cluster D
                clusterFunction <- mainClusterArgs[["clusterFunction"]]
                if (is.character(clusterFunction))
                    clusterFunction <- getBuiltInFunction(clusterFunction)
                
                if(subsample){
                    # Check that mainClustering cluster function takes 
                    # input that is cat
                    # Reason: if subsampling, then the results from subsampling 
                    # sent to the clusterFunction.
                    # This would be caught below anyway, 
                    # but this is a more informative error
                    mainInputType<-"cat"
                    if (!"cat" %in% inputType(clusterFunction))
                        return(
                            "If choosing subsample=TRUE, the clusterFunction used in the mainClustering step must take input that is a categorical (inputType='cat') matrix."
                        )
                }
                else mainInputType<-inputType
                inputMain <- try(.checkCFInput(
                    inputType=mainInputType,
                    desiredInputType = inputType(clusterFunction)),
                    silent=TRUE)
                if(inherits(inputMain, "try-error")){
                    ##In clusterSingle, can build diss if should.
                    mess<-paste("In mainClustering step,",inputMain)
                    if(allowMakeDiss){
                        ## Reason: if subsample=TRUE, then don't want to calculate dissimilarity matrix or anything. Just want check that it is categorical. 
                       if(!subsample){
                           if(mainInputType=="X" & "diss" %in% inputType(clusterFunction)){
                               doDiss<-TRUE 
                               inputMain<-"diss"
                           }
                           else return(mess) 
                       }
                       else return(mess) 
                    }
                    else{
                        return(mess) 
                    }
                } 
                algType <- algorithmType(clusterFunction)
            }
            else{
                return(
                    "Must provide 'clusterFunction' for the mainClustering step to be able to run (give 'clusterFunction' argument via 'mainClusterArgs')"
                )
            }
            #-----
            # removes certain options based on others
            #-----
            reqArgs <- requiredArgs(clusterFunction)
            # remove 'k' if sequential is given
            if (sequential & length(reqArgs) > 0)
                reqArgs <- reqArgs[-which(reqArgs == "k")]
            # remove 'k' if choose 'findBestK=TRUE'
            if (length(reqArgs) > 0 &
                algorithmType(clusterFunction) == "K" &
                "findBestK" %in% names(mainClusterArgs)) {
                if (mainClusterArgs[["findBestK"]])
                    reqArgs <- reqArgs[-which(reqArgs == "k")]
            }
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
                    #----------
                    # check clusterArgs alpha and k are valid
                    #----------
                    if("clusterArgs" %in% names(mainClusterArgs)){
                        clArgs<-mainClusterArgs[["clusterArgs"]]
                        if("alpha" %in% names(clArgs) & algType=="01"){
                            if(is.na(clArgs[["alpha"]]) ||
                                clArgs[["alpha"]]<0 || 
                                clArgs[["alpha"]]>1){
                                return("alpha parameter for type 01 clustering functions must be in (0,1)")
                            }
                        }
                        if("k" %in% names(clArgs) & algType=="K"){
                            if(clArgs[["k"]]<0 || 
                                !is.whole(clArgs[["k"]])){
                                return("k parameter for type K clustering functions must be positive whole number")
                            }
                        }
                    }
            }
            
            
            #--------
            #Check that minSize is valid
            #--------
            if("minSize" %in% names(mainClusterArgs)){
                if (!is.numeric(mainClusterArgs[["minSize"]]) ||
                    mainClusterArgs[["minSize"]] < 0)
                    return(
                        "Invalid value for the 'minSize' parameter in determining the minimum number of samples required in a cluster (in main clustering step)."
                    )
                else
                    mainClusterArgs[["minSize"]] <-
                        round(mainClusterArgs[["minSize"]]) #in case not integer.                
            }
            
            #-------
            ## Check post-processing arguments
            #-------
            doKPostProcess <- FALSE
            
            # whPostProcess args are which of given minClusterArgs are user options to give to .postProcess
            # They will ultimately be the only ones returned back to mainClustering to pass along.
            whPostProcessArgs <-
                which(names(mainClusterArgs) %in% getPostProcessingArgs(clusterFunction))
            if (length(whPostProcessArgs) > 0 | !is.null(mainClusterArgs[["extraArguments"]])) {
                #get rid of wrong args passed because of user confusion between the two
                
                postProcessArgs <- mainClusterArgs[whPostProcessArgs]
                if (clusterFunction@algorithmType == "K") {
                    if ("findBestK"  %in% names(postProcessArgs) &&
                        postProcessArgs[["findBestK"]])
                        doKPostProcess <- TRUE
                    if ("removeSil"  %in% names(postProcessArgs) &&
                        postProcessArgs[["removeSil"]])
                        doKPostProcess <- TRUE
                    if(doKPostProcess){
                        if(inputMain !="diss"){
                            if(!"diss" %in% names(mainClusterArgs)){
                                if(allowMakeDiss & mainInputType=="X"){
                                    doDissPost<-TRUE
                                }
                                else{
                                    return("Cannot do requested post processing (e.g. from arguments 'findBestK' or 'removeSil') without calculating distance matrix")
                                }
                            }
                        }
                    }
                }
                
                
                # remove existing ones, and only keep usable ones.
                # Only relevant in the check run in mainClustering, 
                # where need get rid of extra ones.
                if(!is.null(mainClusterArgs[["extraArguments"]])){
                    exArgsNames<-mainClusterArgs[["extraArguments"]]
                    
                    if (length(whPostProcessArgs) !=
                         length(exArgsNames) & warn){
                             
                             whNotMatch<- !exArgsNames %in%
                                  names(mainClusterArgs)[whPostProcessArgs]
                             warning(
                                 sprintf("Some arguments passed via mainClusterArgs in mainClustering step do not match the algorithmType of the given ClusterFunction object: %s",paste(exArgsNames[whNotMatch],collapse=",")))
                         }
                        
                    whRm<-which(names(mainClusterArgs) %in% exArgsNames)
                    if(length(whRm)>0){
                        mainClusterArgs<-mainClusterArgs[-whRm]
                        mainClusterArgs<-c(mainClusterArgs,postProcessArgs)           
                    }
                    mainClusterArgs[["extraArguments"]] <- names(postProcessArgs)     
                }
            }
            mainClusterArgs[["doKPostProcess"]] <- doKPostProcess
            


            mainClusterArgs[["inputType"]]<-inputMain
			if(!"checkArgs" %in% names(mainClusterArgs[["clusterArgs"]]))
                mainClusterArgs[["clusterArgs"]][["checkArgs"]]<-warn
            mainClusterArgs[["warnings"]]<-warn
            #----------
            # Adapt MainClusteringArgs if using Sequential
            #----------
            if (sequential) {
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
                    #make it exist... not sure if I need this.
                    mainClusterArgs[["clusterArgs"]] <- NULL 
                }
                if(!subsample){
                    # Reason: if subsample=FALSE, and sequential=TRUE 
                    # then need to adjust K in the mainClustering step, so need algorithm of type K
                    if (algorithmType(clusterFunction) != "K") {
                        return(
                            "if subsample=FALSE, sequentical clustering can only be implemented with a clusterFunction for mainClustering step with algorithmType 'K'. See documentation of seqCluster."
                        )
                        #Not sure why this was here!
                        # subsampleArgs <-
                        #                         subsampleArgs[-which(names(subsampleArgs) == "clusterFunction")]
                    }
                    # Reason: subsample=FALSE can't do sequential clustering 
                    # and findBestK=TRUE because need to remove cluster based on
                    # testing many k and finding stable, and if not doing it over 
                    # subsample, then do it over actual clustering
                    if ("findBestK" %in% names(mainClusterArgs)) {
                        if (mainClusterArgs[["findBestK"]])
                            return(
                                "Cannot do sequential clustering where subsample=FALSE and 'findBestK=TRUE' is passed to the mainClustering step via mainClusterArgs. See help documentation of seqCluster."
                            )
                    }     
                }  
            }
        
        }
        
        

        
        #########
        # Checks related to subsample=TRUE
        #########
        if (subsample) {
            if(is.null(subsampleArgs)){
                if(main) subsampleArgs<-list()
                else{
                    #I don't think this can ever actually happen:
                    return("if subsample=TRUE and mainClusterArgs not given, must give argument subsampleArgs so as to identify the clusterFunction")}
            }
            
            # Note if subsample=TRUE, then inputType refers to the subsampling, not the main step.
            inputSub<-inputType
            # Set default for when clusterFunction for subsampling 
            # not given or needs to be overridden (only used if main=TRUE)
            default <- switch(inputSub, 
                              X = "kmeans",
                              cat = "hier01",
                              diss = "pam")
            #---
            # If no cluster function given for subsampling...
            # set some defaults if at top level 
            #---
            if(main & !"clusterFunction" %in% names(subsampleArgs)){
                cfInp<-inputType(mainClusterArgs[["clusterFunction"]])
                if (inputSub %in% cfInp) {
                    mess <-
                        "a clusterFunction was not given for subsampleClustering -- set to be the same as the mainClustering step"
                    if (is.character(mainClusterArgs[["clusterFunction"]]))
                        mess <- sprintf("%s (%s)", mess, mainClusterArgs[["clusterFunction"]])
                    if (warn)
                        .mynote(mess)
                    subsampleArgs[["clusterFunction"]] <- clusterFunction
                    diffSubsampleCF <- FALSE
                }
                else{
                    if(doDiss & inputSub=="X" & "diss" %in% cfInp){
                        # If already going to be calculating 
                        # the dissimilarity matrix for the main step
                        # go ahead and reuse for the subsample step. 
                        mess<-sprintf("a clusterFunction was not given for subsampleClustering -- set to be the same as the mainClustering step")
                        if(warn) .mynote(mess)
                        subsampleArgs[["clusterFunction"]] <- clusterFunction
                        diffSubsampleCF <- FALSE
                    }
                    else{
                        if(warn){
                            mess<-sprintf("a clusterFunction was not given for subsampleClustering (and inputMatrix is not of type allowed by clusterFunction of mainClustering step). The clusterFunction for subsampling was set to the default of %s",
                                default)
                            .mynote(mess)
                        }
                        subsampleArgs[["clusterFunction"]] <- default
                        diffSubsampleCF <- TRUE
                    }
                }
            }
            
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
                inputSub<-try(.checkCFInput(inputType=inputSub, 
                    desiredInputType=inputType(subsampleCF)),
                    silent=TRUE)
                if(inherits(inputSub, "try-error")){
                    mess<-paste("In subsampling clustering step,",inputSub)
                    ##In clusterSingle, can build diss if should.
                    if(allowMakeDiss){
                       if(inputType=="X" & "diss" %in% inputType(subsampleCF)){
                           doDiss<-TRUE 
                           inputSub<-"diss"
                       }
                       else return(mess) 
                    }
                    else{
                        return(mess) 
                    }
                } 
                
                # Check that subsample clustering function is of type 'k' if sequential=TRUE
                # Reason: seqCluster requires subsampling cluster function to be of type "K"
                if (sequential & algorithmType(subsampleCF) != "K") {
                    if (warn)
                        warning(
                            sprintf(
                                "If subsample=TRUE, sequentical clustering can only be implemented with a clusterFunction for subsampling that has algorithmType 'K'. See documentation of seqCluster. Will ignore 'clusterFunction' argument of subsampleArgs and set to default of %s",
                                default
                            )
                        )
                    subsampleArgs[["clusterFunction"]] <- default
                    diffSubsampleCF <- TRUE
                }
            }
            else return("Must give 'clusterFunction' value to 'subsampleClustering'.")
            #-----
            #Check the classifyMethod argument
            #-----
            if(!is.null(subsampleArgs[["classifyMethod"]])){
                if (is.null(subsampleCF@classifyFUN)) {
                    if ("classifyMethod" %in% names(subsampleArgs) &&
                        subsampleArgs[["classifyMethod"]] != "InSample")
                        if(warn) warning(
                            "Cannot set 'classifyMethod' to anything but 'InSample' if do not specify a clusterFunction in subsampleArgs that has a non-null classifyFUN slot. Changing classifyMethod to 'InSample'."
                        )
                    subsampleArgs[["classifyMethod"]] <- "InSample"
                }
                if(subsampleArgs[["classifyMethod"]]!="InSample"){
                    if(!inputSub %in% subsampleCF@inputClassifyType){
                        if(warn) warning(sprintf("The clusterFunction for subsampling chosen does not provide ability to classify new samples from input type %s to clusters in the subsampling step (needed by classifyMethod=%s in 'subsampleArgs'). Changing classifyMethod to 'InSample'.", inputSub, subsampleArgs[["classifyMethod"]]))
                        subsampleArgs[["classifyMethod"]]<-"InSample"
                    }
                }
            }

            
            ##------
            ## Check have required args for subsample. 
            ## If missing, 'borrow' those args from mainClusterArgs.
            ##------
            reqSubArgs <- requiredArgs(subsampleCF)
            
            #Reason: sequential sets k for the subsampling via k0
            if (sequential & length(reqSubArgs) > 0)
                reqSubArgs <- reqSubArgs[-which(reqSubArgs == "k")]
            if (length(reqSubArgs) > 0){
                if(main){
                    #check if can borrow...
                    if ("clusterArgs" %in% names(mainClusterArgs)) {
                        mainReqArgs <- requiredArgs(clusterFunction)
                        mainReqArgs <-
                            mainReqArgs[mainReqArgs %in%
                            names(mainClusterArgs[["clusterArgs"]])]
                        if (!is.null(subsampleArgs) &&
                            "clusterArgs" %in% names(subsampleArgs)) {
                            # check if existing clusterArgs has required names already
                            # Define missingArgs to be those that are needed and available from main 
                            if (!all(reqSubArgs %in%
                                names(subsampleArgs[["clusterArgs"]]))) {
                                missingArgs <-
                                    reqSubArgs[!reqSubArgs %in%
                                    names(subsampleArgs[["clusterArgs"]])]
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
            
            #---
            #Check don't give distance function and subsample=TRUE
            #Reason, if subsample=TRUE, user can't set distance function because use diss provided from subsampling.
            #---
            if(main &&
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
            subsampleArgs[["inputType"]]<-inputSub
            if(!"checkArgs" %in% names(subsampleArgs[["clusterArgs"]]))
                subsampleArgs[["clusterArgs"]][["checkArgs"]]<-warn
            subsampleArgs[["warnings"]]<-warn
						
        }    
        
        #########
        # Check seqArgs
        #########
        if(sequential){
            if(is.null(seqArgs))
                return("if sequential=TRUE, must give seqArgs so as to identify k0 and beta")
            
            ## Test parameters that need to be whole positive numbers
            testFUN<-function(name){
                if(name %in% seqArgs){
                    if(is.na(seqArgs[[name]]) || !is.whole(seqArgs[[name]]) || seqArgs[[name]]<0) return(sprintf("%s must be a positive whole number",name))
                    else return(TRUE)                    
                }
                else return(TRUE)
            }
            if(!"k0"%in%names(seqArgs)) {
                return("required argument 'k0' is missing for the sequential clustering step")
            }
            else{
                out<-testFUN("k0")
                if(!is.logical(out)) return(out)
            }
            if(!"beta" %in% names(seqArgs))
                return("required argument 'beta' is missing for the sequential clustering step")
            else{
                if(is.na(seqArgs[["beta"]]) || seqArgs[["beta"]]<=0 || seqArgs[["beta"]]>=1 ) return("beta parameter for sequential clustering step must be in (0,1)")
            }
            for(nam in c("top.can","remain.n","k.min","k.max")){
                out<-testFUN(nam)
                if(!is.logical(out)) return(out)
            }
            seqArgs[["warnings"]]<-warn
        }    
        return(
            list(
                mainClusterArgs = mainClusterArgs,
                subsampleArgs = subsampleArgs,
                seqArgs = seqArgs,
                doDiss=doDiss,
                doDissPost=doDissPost
            )
        )
        
    }


# from Martin Maechler  to test if whole number https://stat.ethz.ch/pipermail/r-help/2003-April/032471.html
is.whole <- function(a, tol = 1e-7) { 
   is.eq <- function(x,y) { 
	 r <- all.equal(x,y, tol=tol)
	 is.logical(r) && r 
   }
   (is.numeric(a) && is.eq(a, floor(a))) ||
   (is.complex(a) && {ri <- c(Re(a),Im(a)); is.eq(ri, floor(ri))})
}