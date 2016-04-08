#' Create a matrix of clustering across values of k
#'
#' Given a range of k's, this funciton will return a matrix with the clustering of the samples
#' across the range, which can be passed to \code{plotTracking} for visualization.
#'
#' @param data the data on which to run the clustering. Must either be a data.frame/matrix (with samples in rows) or a list of datasets overwhich the clusterings should be run.
#' @param ks the range of k values (see details for meaning for different choices). 
#' @param alphas values of alpha to be tried. Only used for subsampleClusterMethod either 'tight' or 'hierarchical'.
#' @param findBestK values of findBestK to be tried (logical) (only for 'pam').
#' @param sequential values of sequential to be tried (logical) (only for 'pam')
#' @param removeSil values of removeSil to be tried (logical) (only for 'pam')
#' @param subsample values of subsample to be tried (logical).
#' @param silCutoff values of silCutoff to be tried (only for 'pam')
#' @param clusterMethod method used in clustering of subsampled data passed to argument 'cluserFunction' of \code{\link{clusterD}}. Note that unlike other functions of this package, this must be a character vector of pre-defined clustering techniques provided by the package, and can not be user-defined.
#' @param clusterDArgs list of arguments to be passed to \code{\link{clusterD}}
#' @param subsampleArgs list of arguments to be passed to \code{\link{subsampleClustering}}
#' @param seqArgs list of arguments to be passed to \code{\link{seqCluster}}
#' @param ncores the number of threads
#' @param random.seed a value to set seed before each run of clusterAll (so that all of the runs are run on the same subsample of the data)
#' @param run logical. If FALSE, doesn't run clustering, but just returns matrix of parameters that will be run for inspection by user (with rownames equal to the names of the resulting column names of clMat object that would be returned if \code{run=TRUE}). 
#' @param paramMatrix matrix or data.frame. If given, the algorithm will bypass creating the matrix of possible parameters, and will use the given matrix. There are basically no checks as to whether this matrix is in the right format, and is only intended to be used to feed the results of setting \code{run=FALSE} back into the algorithm (see example). 
#' @param ... arguments to be passed on to mclapply (if ncores>1)
#'
#' @details While the function allows for multiple values of clusterMethod, the code does not reuse the same subsampling matrix and try different clusterMethods on it. If sequential=TRUE, different subsampleClusterMethods will create different sets of data to subsample so it is not possible; if sequential=FALSE, we have not implemented functionality for this reuse. Setting the \code{random.seed} value, however, should mean that the subsampled matrix is the same for each, but there is no gain in computational complexity (i.e. each subsampled co-occurence matrix is recalculated for each set of parameters). 
#'
#' @details Note that the behavior of compareChoices for dimensionality reduction is slightly different if the input is a list of datasets rather than a matrix. If the input is a single matrix, a single dimensionality step is performed, while if the input is a list of datasets, the dimensionality reduction step is performed for every combination (i.e. the program is not smart in realizing it is the same set of data across different dimensions and so only one dimensionality reduction is needed). 
#'
#' @details The argument 'ks' is interpreted differently for different choices of the other parameters. When/if sequential=TRUE, ks defines the argument k0 of \code{\link{seqCluster}}. When/if clusterMethod="pam" and "findBestK=TRUE", ks defines the kRange argument of \code{\link{clusterD}} unless kRange is specified by the user via the clusterDArgs; note this means that the default option of setting kRange that depends on the input k (see \code{\link{clusterD}}) is not available in compareChoices. 
#' @return If \code{run=TRUE}, a list with the following objects:
#' \itemize{
#' \item{\code{clMat}}{a matrix of with each row corresponding to a clustering and each column a sample.}
#' \item{\code{clusterInfo}}{a list with information regarding clustering result (only relevant entries for those clusterings with sequential=TRUE)}
#' }
#' @return If \code{run=FALSE}, a matrix which each row corresponding to a set of parameters to pass to \code{\link{clusterAll}}. This matrix can be given to \code{paramMatrix} to then run these parameter choices.
#'
#' @examples
#' data(simData)
#' #clustering using pam: try using different dimensions of pca and different k
#' ps<-c(5,10,50)
#' names(ps)<-paste("npc=",ps,sep="")
#' pcaData<-stats::prcomp(simData, center=TRUE, scale=TRUE)
#' 
#' #check how many and what runs user choices will imply:
#' #Note that this causes error:
#' \dontrun{
#' checkParams <- compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}), clusterMethod="pam",
#' ks=2:4,findBestK=c(TRUE,FALSE),run=FALSE)
#' }
#' #fixes error, but really, not clear what best subsampling k should be
#' checkParams <- compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}), clusterMethod="pam",
#' ks=2:4,findBestK=c(TRUE,FALSE),run=FALSE,subsampleArgs=list("k"=3))
#' #Now actually run it
#' cl <- compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}), 
#' clusterMethod="pam",ks=2:4,findBestK=c(TRUE,FALSE),
#' subsampleArgs=list("k"=3))
#' colnames(cl$clMat) 
#' #make names shorter for plotting
#' colnames(cl$clMat)<-gsub("TRUE","T",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("FALSE","F",colnames(cl$clMat))
#' colnames(cl$clMat)<-gsub("k=NA,","",colnames(cl$clMat))
#' par(mar=c(2,10,1,1))
#' plotTracking(cl$clMat,axisLine=-2)
#' #get rid of some of the choices manually
#' checkParams<-checkParams[-c(1,2),]
#' clSmaller<-compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}),
#' paramMatrix=checkParams)
#' 
#' 
#'
#' \dontrun{
#'	#following code takes around 1+ minutes to run because of the subsampling that is redone each time:
#'	system.time(clusterTrack<-compareChoices(simData, ks=2:15, 
#'	alphas=c(0.1,0.2,0.3), findBestK=c(TRUE,FALSE),sequential=c(FALSE),
#'	subsample=c(FALSE),removeSil=c(TRUE), clusterMethod="pam", 
#'	clusterDArgs = list(minSize = 5,kRange=2:15),ncores=1,random.seed=48120))

#' }
#' 
#Work up example:
# clusterTrack<-compareChoices(simData, ks=2:3,
# alphas=c(0.1), findBestK=c(TRUE),sequential=c(FALSE),
# subsample=c(TRUE),removeSil=c(TRUE), clusterMethod=c("pam","tight","hierarchical",
# clusterDArgs = list(minSize = 5,kRange=2:15),subsampleArgsncores=1,random.seed=48120)

#' @rdname compareChoices
setMethod(
  f = "compareChoices",
  signature = signature(x = "matrix"),
  definition = function(x, ks=3:5, clusterMethod="pam", alphas=0.1, findBestK=FALSE,sequential=FALSE,
                        removeSil=FALSE, subsample=FALSE,silCutoff=0,
                        dimReduce="none",nVarDims=NA,nPCADims=NA,
                        clusterDArgs=list(minSize=5),
                        subsampleArgs=list(resamp.num=50),
                        seqArgs=list(beta=0.9,k.min=3, verbose=FALSE),
                        transFun=NULL,isCount=FALSE,
                        ncores=1,random.seed=NULL,run=TRUE,paramMatrix=NULL,...
  ){
    origX<-x
    transObj<-.transData(x,nPCADims=nPCADims, nVarDims=nVarDims,dimReduce=dimReduce,transFun=transFun,isCount=isCount)
    x<-transObj$x
    if(!is.null(dim(x)) && NCOL(x)!=NCOL(origX)) stop("Error in the internal transformation of x")
    transFun<-transObj$transFun #need it later to create clusterCellsObject
    
    if(!is.null(dim(x))) x<-list(dataset1=x) #if npcs=NA, then .transData returns a matrix.
    outval<-compareChoices(x, ks=ks,clusterMethod=clusterMethod,alphas=alphas,findBestK=findBestK,
                           sequential=sequential,removeSil=removeSil,subsample=subsample,silCutoff=silCutoff,
                           clusterDArgs=clusterDArgs,subsampleArgs=subsampleArgs,seqArgs=seqArgs,ncores=ncores,
                           #dimReduce="none",nVarDims=NA,nPCADims=NA, #because already did it above
                           random.seed=random.seed,run=run,paramMatrix=paramMatrix,...)
    ##########
    ## Convert to clusterCells Object
    ##########
    if(run){#browser()
      retval <- clusterCells(origX, outval$clMat[,1], transformation=transFun)
      retval@clusterLabels<-outval$clMat
      retval@clusterInfo<-outval$clusterInfo
      retval@clusterType <- rep("compareChoices",NCOL(outval$clMat))
      validObject(retval)
      return(retval)
    }
    else return(outval)
  }

)

#' @rdname compareChoices
setMethod(
  f = "compareChoices",
  signature = signature(x = "list"),
  definition = function(x, ks, clusterMethod, alphas=0.1, findBestK=FALSE,sequential=FALSE,
                        removeSil=FALSE, subsample=FALSE,silCutoff=0,
                        clusterDArgs=list(minSize=5),
                        subsampleArgs=list(resamp.num=50),
                        seqArgs=list(beta=0.9,k.min=3, verbose=FALSE),
                        ncores=1,random.seed=NULL,run=TRUE,paramMatrix=NULL,...
  )
  {
    data<-x
    if(!all(sapply(data,function(y){is.matrix(y) || is.data.frame(y)}))) stop("if data is a list, it must be a list with each element of the list a data.frame or matrix")
    #check all same number of observations:
    if(!length(unique(sapply(data,NCOL)))==1) stop("All data sets must have the same number of observations")
    if(is.null(names(data))) names(data)<-paste("dataset",1:length(data),sep="")
    dataList<-data
    dataName<-names(dataList)
    if(is.null(paramMatrix)){
      
      # #check if ks=NA; no longer an option
      # if(is.na(ks)){
      # 	if(!all(clusterMethod=="pam")) stop("if clusterMethod includes methods that are not 'pam', must provide ks")
      # 	else if(!all(findBestK)) stop("if findBestK!=TRUE, must provide ks")
      # }
      param<-expand.grid(dataset=dataName,#dimReduce="none",nVarDims=NA,nPCADims=NA,
                         k=ks,alpha=alphas,findBestK=findBestK,sequential=sequential,
                         removeSil=removeSil,subsample=subsample,
                         clusterMethod=clusterMethod,silCutoff=silCutoff)
      #don't vary them across ones that don't matter (i.e. 0-1 versus K); 
      #code sets to single value and then will do unique
      #also deals with just in case the user gave duplicated values of something by mistake.
      typeK<-which(param[,"clusterMethod"] %in% c("pam"))
      if(length(typeK)>0){
        param[typeK,"alpha"]<-NA #just a nothing value
        whFindBestK<-which(param[,"findBestK"])
        if(length(whFindBestK)>0){ #remove 'k' and see if same
          if(!"kRange" %in% names(clusterDArgs)) clusterDArgs[["kRange"]]<-ks
          #if findBestK=TRUE, and sequential=FALSE, then user needs to set k via subsampleArgs
          whNoSeq<-which(!param[,"sequential"])
          if(length(intersect(whFindBestK,whNoSeq))>0){
            param[intersect(whFindBestK,whNoSeq),"k"]<-NA
            if(is.null(subsampleArgs[["k"]])) stop("must provide k in subsampleArgs for those with findBestK=TRUE and sequential=FALSE")
            #        else param[intersect(whFindBestK,whNoSeq),"k"]<-subsampleArgs[["k"]]
          } 
        }
      }
      type01<-which(param[,"clusterMethod"] %in% c("hierarchical","tight"))
      if(length(type01)>0){
        param[type01,"findBestK"]<-FALSE
        param[type01,"removeSil"]<-FALSE
        param[type01,"silCutoff"]<-0
      }

      #if provide kRange in clusterDArgs, and findBestK=TRUE, don't need to search over different k
      param<-unique(param)  

      #deal with those that are invalid combinations:
      whInvalid<-which(!param[,"subsample"] & param[,"sequential"] & param[,"findBestK"])
      if(length(whInvalid)>0) param<-param[-whInvalid,]
      
      whExtra<-which(!param[,"subsample"] & param[,"sequential"] & param[,"findBestK"])
      if(length(whInvalid)>0) param<-param[-whInvalid,]
      
      if(nrow(param)<=1) stop("set of parameters imply only 1 combination")
      
      #find names of the parameter combinations.
      whVary<-which(apply(param,2,function(x){length(unique(x))>1}))
      if(length(whVary)>0) cnames<-apply(param[,whVary,drop=FALSE],1,function(x){
        paste(colnames(param)[whVary],x,sep="=",collapse=",")})
      else stop("set of parameters imply only 1 combination")
      
      cnames<-gsub("dataset=","",cnames)
      cnames<-gsub("= ","=",cnames)
      cnames[param[,"sequential"]]<-gsub("k=","k0=",cnames[param[,"sequential"]])
      cat(nrow(param),"parameter combinations,",sum(param[,"sequential"]),"use sequential method.\n")
      rownames(param)<-cnames
    }
    else{
      if(!run) stop("If paramMatrix is given, run should be TRUE. Otherwise there is no effect.")
      if(is.null(paramMatrix)) stop("invalid input for paramMatrix; must be data.frame or matrix")
      param<-paramMatrix
      if(is.null(rownames(paramMatrix))) stop("input paramMatrix must have row names")
      cnames<-rownames(paramMatrix)

    }
    

    paramFun<-function(i){
      par<-param[i,]
      #make them logical values... otherwise adds a space before the TRUE and doesn't recognize.
      #well, sometimes. Maybe don't need this?
      removeSil<-as.logical(gsub(" ","",par["removeSil"]))
      sequential<-as.logical(gsub(" ","",par["sequential"]))
      subsample<-as.logical(gsub(" ","",par["subsample"]))
      findBestK<-as.logical(gsub(" ","",par["findBestK"]))
      clusterMethod<-as.character(par[["clusterMethod"]])
      if(!is.na(par[["k"]])){
        if(sequential) seqArgs[["k0"]]<-par[["k"]] 
        else{
          #to be safe, set both in case user set one. 
          subsampleArgs[["k"]]<-par[["k"]]
          clusterDArgs[["k"]]<-par[["k"]]
        }			
      }
      clusterDArgs[["alpha"]]<-par[["alpha"]]
      clusterDArgs[["findBestK"]]<-findBestK
      clusterDArgs[["removeSil"]]<-removeSil
      clusterDArgs[["silCutoff"]]<-par[["silCutoff"]]
      clusterDArgs[["checkArgs"]]<-FALSE #turn off printing of warnings that arguments off
      if(!is.null(random.seed)) set.seed(random.seed)
      clusterAll(x=dataList[[par[["dataset"]]]],  subsample=subsample,clusterFunction=clusterMethod,  clusterDArgs=clusterDArgs,subsampleArgs=subsampleArgs,
                 seqArgs=seqArgs, sequential=sequential,transFun=function(x){x}) #dimReduce=dimReduce,ndims=ndims,
    }
    if(run){
      cat("Running Clustering on Parameter Combinations...")
      if(ncores>1){
        out<-mclapply(1:nrow(param),FUN=paramFun,mc.cores=ncores,...)
        nErrors<-which(sapply(out,function(x){inherits(x, "try-error")}))
        if(length(nErrors)>0)stop(nErrors,"parameter values hit an error. The first was:\n",out[nErrors[1]])
      }
      else out<-lapply(1:nrow(param),FUN=paramFun)
      cat("done.\n")
      clMat<-sapply(out,function(x){primaryCluster(x)})
      
      colnames(clMat)<-unname(cnames)
      pList<-lapply(1:nrow(param),function(i){
        x<-param[i,]
        names(x)<-colnames(param)
        return(x)})
      clInfo<-mapply(pList,out,FUN=function(x,y){
        c(list(choicesParam=x),clusterInfo(y))
      },SIMPLIFY=FALSE)
      
      return(list(clMat=clMat,clusterInfo=clInfo,paramMatrix=param))
    }
    else{
      cat("Returning Parameter Combinations without running them (to run them choose run=TRUE)\n")
      return(param)
    }
  }
)

#Not yet implemented
#' @rdname compareChoices
setMethod(
  f = "compareChoices",
  signature = signature(x = "ClusterCells"),
  definition = function(x, ks=3:5, clusterMethod, alphas=0.1, findBestK=FALSE,sequential=FALSE,
                        removeSil=FALSE, subsample=FALSE,silCutoff=0,
                        dimReduce="none",nVarDims=NA,nPCADims=NA,
                        clusterDArgs=list(minSize=5),
                        subsampleArgs=list(resamp.num=50),
                        seqArgs=list(beta=0.9,k.min=3, verbose=FALSE),
                        ncores=1,random.seed=NULL,run=TRUE,paramMatrix=NULL,...
  )
  {
    
    outval<-compareChoices(assay(x), ks=ks, clusterMethod=clusterMethod, alphas=alphas,
                           findBestK=findBestK,sequential=sequential,
                           removeSil=removeSil, subsample=subsample,silCutoff=silCutoff,
                           dimReduce=dimReduce,nVarDims=nVarDims,nPCADims=nPCADims,
                           clusterDArgs=clusterDArgs,
                           subsampleArgs=subsampleArgs,
                           seqArgs=seqArgs,
                           transFun=transformation(x),
                           ncores=ncores,random.seed=random.seed,run=run,paramMatrix=paramMatrix,eraseOld=FALSE,...
    )
    if(run){
      ##Check if compareChoices already ran previously
      ppIndex<-pipelineClusterIndex(x,print=FALSE)
      if(!is.null(ppIndex)){ #need to change the clusterType values (or erase them) before get new ones
        if(eraseOld){ #remove all of them, not just current
          x<-removeClusters(x,ppIndex[,"index"])
        }
        else{
          if(0 %in% ppIndex[,"iteration"]){
            newIteration<-max(ppIndex[,"iteration"])+1
            whCurrent<-ppIndex[ppIndex[,"iteration"]==0,"index"]
            updateCluster<-clusterType(x)
            updateCluster[whCurrent]<-paste(updateCluster[whCurrent],newIteration,sep="_")
            x@clusterType<-updateCluster          
          }
        }
        
      }
      
      retval<-addClusters(outval,x)
      validObject(retval)
      return(retval)
    }
    else return(outval)
  }
)


#' @rdname compareChoices
setMethod(
  f = "compareChoices",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, ks=3:5, clusterMethod, alphas=0.1, findBestK=FALSE,sequential=FALSE,
                        removeSil=FALSE, subsample=FALSE,silCutoff=0,
                        dimReduce="none",nVarDims=NA,nPCADims=NA,
                        clusterDArgs=list(minSize=5),
                        subsampleArgs=list(resamp.num=50),
                        seqArgs=list(beta=0.9,k.min=3, verbose=FALSE),
                        transFun=NULL,isCount=FALSE,
                        ncores=1,random.seed=NULL,run=TRUE,paramMatrix=NULL,...
  )
  {
    outval<-compareChoices(assay(x), ks=ks, clusterMethod=clusterMethod, alphas=alphas,
             findBestK=findBestK,sequential=sequential,
             removeSil=removeSil, subsample=subsample,silCutoff=silCutoff,
             dimReduce=dimReduce,nVarDims=nVarDims,nPCADims=nPCADims,
             clusterDArgs=clusterDArgs,
             subsampleArgs=subsampleArgs,
             seqArgs=seqArgs,
             transFun=transFun,isCount=isCount,
             ncores=ncores,random.seed=random.seed,run=run,paramMatrix=paramMatrix,...
    )
    if(run){
      retval <- clusterCells(x, primaryCluster(outval), transformation(outval))
      retval@clusterLabels<-outval@clusterLabels
      retval@clusterInfo <- clusterInfo(outval)
      retval@clusterType <- clusterType(outval) 
      validObject(retval)
      return(retval)
      
    }  #need to redo it to make sure get any other part of summarized experiment
    else return(outval)  
  }
)



