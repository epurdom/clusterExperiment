#' Create a matrix of clustering across values of k
#'
#' Given a range of k's, this funciton will return a matrix with the clustering of the samples
#' across the range, which can be passed to \code{plotTracking} for visualization.
#'
#' @param x the data on which to run the clustering. Can be: data.frame/matrix (with samples in rows), a list of datasets overwhich the clusterings should be run, a \code{SummarizedExperiment} object, or a \code{ClusterCells} object.
#' @param ks the range of k values (see details for meaning for different choices). 
#' @param alphas values of alpha to be tried. Only used for subsampleClusterMethod either 'tight' or 'hierarchical'.
#' @param findBestK values of findBestK to be tried (logical) (only for 'pam').
#' @param sequential values of sequential to be tried (logical) (only for 'pam')
#' @param removeSil values of removeSil to be tried (logical) (only for 'pam')
#' @param subsample values of subsample to be tried (logical).
#' @param silCutoff values of silCutoff to be tried (only for 'pam')
#' @param clusterMethod method used in clustering of subsampled data passed to argument 'cluserFunction' of \code{\link{clusterD}}. Note that unlike other functions of this package, this must be a character vector of pre-defined clustering techniques provided by the package, and can not be user-defined.
#' @param dimReduce character vector of dimensionality reduction methods to try, that should be some combination of 'PCA', 'mostVar' and 'none'
#' @param nVarDims vector of the number of the most variable features to keep (when "mostVar" is identified in \code{dimReduce}). If NA is included, then the full dataset will also be included.
#' @param nPCADims vector of the number of PCs to use (when 'PCA' is identified in \code{dimReduce}). If NA is included, then the full dataset will also be included.
#' @param transFun function to use to transform the data
#' @param isCount logical. If transFun is missing, will be used to determine the transformation. log(x+1) will be the transformation for isCount=TRUE and otherwise the identify function x.
#' @param eraseOld logical. Only relevant if input \code{x} is of class \code{ClusterCells}. If TRUE, will erase existing pipeline results (compareChoices as well as mergeClusters and findSharedClusters). If FALSE, existing pipeline results will have "\code{_i}" added to the clusterType value, where \code{i} is one more than the largest such existing pipeline clusterType.
#' @param clusterDArgs list of arguments to be passed to \code{\link{clusterD}}
#' @param subsampleArgs list of arguments to be passed to \code{\link{subsampleClustering}}
#' @param seqArgs list of arguments to be passed to \code{\link{seqCluster}}
#' @param ncores the number of threads
#' @param random.seed a value to set seed before each run of clusterAll (so that all of the runs are run on the same subsample of the data)
#' @param run logical. If FALSE, doesn't run clustering, but just returns matrix of parameters that will be run, for the purpose of inspection by user (with rownames equal to the names of the resulting column names of clMat object that would be returned if \code{run=TRUE}). 
#' @param paramMatrix matrix or data.frame. If given, the algorithm will bypass creating the matrix of possible parameters, and will use the given matrix. There are basically no checks as to whether this matrix is in the right format, and is only intended to be used to feed the results of setting \code{run=FALSE} back into the algorithm (see example). 
#' @param ... arguments to be passed on to mclapply (if ncores>1)
#'
#' @details While the function allows for multiple values of clusterMethod, the code does not reuse the same subsampling matrix and try different clusterMethods on it. If sequential=TRUE, different subsampleClusterMethods will create different sets of data to subsample so it is not possible; if sequential=FALSE, we have not implemented functionality for this reuse. Setting the \code{random.seed} value, however, should mean that the subsampled matrix is the same for each, but there is no gain in computational complexity (i.e. each subsampled co-occurence matrix is recalculated for each set of parameters). 
#'
#' @details Note that the behavior of compareChoices for dimensionality reduction is slightly different if the input is a list of datasets rather than a matrix. If the input is a single matrix, a single dimensionality step is performed, while if the input is a list of datasets, the dimensionality reduction step is performed for every combination (i.e. the program is not smart in realizing it is the same set of data across different dimensions and so only one dimensionality reduction is needed). 
#'
#' @details The argument 'ks' is interpreted differently for different choices of the other parameters. When/if sequential=TRUE, ks defines the argument k0 of \code{\link{seqCluster}}. When/if clusterMethod="pam" and "findBestK=TRUE", ks defines the kRange argument of \code{\link{clusterD}} unless kRange is specified by the user via the clusterDArgs; note this means that the default option of setting kRange that depends on the input k (see \code{\link{clusterD}}) is not available in compareChoices. 
#' @return If \code{run=TRUE} and the input either a matrix, a \code{SummarizedExperiment} object, or a \code{ClusterCells} object, will return a \code{ClusterCells} Object, where the results are stored as clusterings with clusterType \code{compareChoices}. Depending on \code{eraseOld} argument above, this will either delete existing such objects, or change the clusterType of existing objects. See argument \code{eraseOld} above. Arbitrarily the first clustering is set as the primaryClusteringIndex
#' 
#' @return If \code{run=TRUE} and the input is a list of data sets, a list with the following objects:
#' \itemize{
#' \item{\code{clMat}}{a matrix of with each row corresponding to a clustering and each column a sample.}
#' \item{\code{clusterInfo}}{a list with information regarding clustering result (only relevant entries for those clusterings with sequential=TRUE)}
#' \item{\code{paramMatrix}}{a matrix giving the parameters of each clustering, where each column is a possible parameter set by the user and passed to \code{\link{clusterAll}} and and each row of paramMatrix corresponds to a clustering in \code{clMat} }
#' \item{\code{clusterDArgs}}{a list of (possibly modified) arguments to clusterDArgs}
#' \item{\code{seqArgs=seqArgs}}{a list of (possibly modified) arguments to seqArgs}
#' \item{\code{subsampleArgs}}{a list of (possibly modified) arguments to subsampleArgs}
#' }
#' @return If \code{run=FALSE} a list similar to that described above, but without the clustering results.
#'
#' @examples
#' data(simData)
#' #Example: clustering using pam with different dimensions of pca and different k and whether remove negative silhouette values

#' #check how many and what runs user choices will imply:
#' checkParams <- compareChoices(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterMethod="pam",
#' ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE),run=FALSE)
#' print(checkParams$paramMatrix)
#' #Now actually run it
#' cl <- compareChoices(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterMethod="pam",ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE))
#' print(cl)
#' colnames(allClusters(cl)) 
#' #make names shorter for plotting
#' clMat<-allClusters(cl)
#' colnames(clMat)<-gsub("TRUE","T",colnames(clMat))
#' colnames(clMat)<-gsub("FALSE","F",colnames(clMat))
#' colnames(clMat)<-gsub("k=NA,","",colnames(clMat))
#' par(mar=c(2,10,1,1))
#' plotTracking(clMat,axisLine=-2)
#' #get rid of some of the choices manually
#' #note that the supplement arguments could have been changed too, so
#' #we give those to compareChoices as well.
#' checkParamsMat<-checkParams$paramMatrix[-c(1,2),]
#' clSmaller<-compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}),
#' paramMatrix=checkParamsMat,subsampleArgs=checkParam$subsampleArgs,
#' seqArgs=checkParam$seqArgs,clusterDArgs=checkParam$clusterDArgs)
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
  definition = function(x, 
                        dimReduce="none",nVarDims=NA,nPCADims=NA,
                        transFun=NULL,isCount=FALSE,
                        ...
  ){
    origX<-x
    transObj<-.transData(x,nPCADims=nPCADims, nVarDims=nVarDims,dimReduce=dimReduce,transFun=transFun,isCount=isCount)
    x<-transObj$x
    if(!is.null(dim(x)) && NCOL(x)!=NCOL(origX)) stop("Error in the internal transformation of x")
    transFun<-transObj$transFun #need it later to create clusterCellsObject
    
    if(!is.null(dim(x))) x<-list(dataset1=x) #if npcs=NA, then .transData returns a matrix.
    outval<-compareChoices(x, ...)
    ##########
    ## Convert to clusterCells Object
    ##########
    if("clMat" %in% names(outval)){
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
      param<-expand.grid(dataset=dataName,#dimReduce="none",nVarDims=NA,nPCADims=NA,
                         k=ks,alpha=alphas,findBestK=findBestK,sequential=sequential,
                         removeSil=removeSil,subsample=subsample,
                         clusterMethod=clusterMethod,silCutoff=silCutoff)
      ###########
      #Check param matrix:
      #don't vary them across ones that don't matter (i.e. 0-1 versus K); 
      #code sets to single value and then will do unique
      #also deals with just in case the user gave duplicated values of something by mistake.
      ###########
      typeK<-which(param[,"clusterMethod"] %in% c("pam"))
      if(length(typeK)>0){
        param[typeK,"alpha"]<-NA #just a nothing value, because doesn't mean anything here
        
        #if findBestK make sure other arguments make sense:
        whFindBestK<-which(param[,"findBestK"])
        if(length(whFindBestK)>0){ 
          #by default make kRange in clusterD equal to the ks. Note this will be true of ALL
          if(!"kRange" %in% names(clusterDArgs)) clusterDArgs[["kRange"]]<-ks 
          
          #if findBestK=TRUE, and sequential=FALSE, then need to set 'k'=NA
          whNoSeq<-which(!param[,"sequential"])
          if(length(intersect(whFindBestK,whNoSeq))>0){
            param[intersect(whFindBestK,whNoSeq),"k"]<-NA
          }
          
          #and if subsample=TRUE, then user needs to set k via subsampleArgs
          whNoSeqSub<-which(!param[,"sequential"] & param[,"subsample"])
          if(length(intersect(whFindBestK,whNoSeqSub))>0 & is.null(subsampleArgs[["k"]])) stop("must provide k in subsampleArgs because there are combinations of findBestK=TRUE, sequential=FALSE and subsample=TRUE. (Note this will set 'k' for all that subsample, even for other parameter combinations)")
        }
      }
      type01<-which(param[,"clusterMethod"] %in% c("hierarchical","tight"))
      if(length(type01)>0){
        param[type01,"findBestK"]<-FALSE
        param[type01,"removeSil"]<-FALSE
        param[type01,"silCutoff"]<-0
      }
      param<-unique(param)  

      #####
      #deal with those that are invalid combinations:
      #####
      whInvalid<-which(!param[,"subsample"] & param[,"sequential"] & param[,"findBestK"])
      if(length(whInvalid)>0) param<-param[-whInvalid,]
      
      whExtra<-which(!param[,"subsample"] & param[,"sequential"] & param[,"findBestK"])
      if(length(whInvalid)>0) param<-param[-whInvalid,]
      
      if(nrow(param)<=1) stop("set of parameters imply only 1 combination")
      
      #give names to the parameter combinations.
      whVary<-which(apply(param,2,function(x){length(unique(x))>1}))
      if(length(whVary)>0) cnames<-apply(param[,whVary,drop=FALSE],1,function(x){
        paste(colnames(param)[whVary],x,sep="=",collapse=",")})
      else stop("set of parameters imply only 1 combination")
      
      cnames<-gsub("dataset=","",cnames)
      cnames<-gsub("= ","=",cnames)
      cnames[param[,"sequential"]]<-gsub("k=","k0=",cnames[param[,"sequential"]])
      rownames(param)<-cnames
    }
    else{
      if(!run) stop("If paramMatrix is given, run should be TRUE. Otherwise there is no effect.")
      if(is.null(paramMatrix)) stop("invalid input for paramMatrix; must be data.frame or matrix")
      param<-paramMatrix
      if(is.null(rownames(paramMatrix))) stop("input paramMatrix must have row names")
      cnames<-rownames(paramMatrix)

    }
    
    cat(nrow(param),"parameter combinations,",sum(param[,"sequential"]),"use sequential method.\n")
      
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
      
      return(list(clMat=clMat,clusterInfo=clInfo,paramMatrix=param,clusterDArgs=clusterDArgs,seqArgs=seqArgs,subsampleArgs=subsampleArgs))
    }
    else{
       cat("Returning Parameter Combinations without running them (to run them choose run=TRUE)\n")
      return(list(paramMatrix=param,clusterDArgs=clusterDArgs,seqArgs=seqArgs,subsampleArgs=subsampleArgs))
    }
  }
)

#' @rdname compareChoices
setMethod(
  f = "compareChoices",
  signature = signature(x = "ClusterCells"),
  definition = function(x, dimReduce="none",nVarDims=NA,nPCADims=NA,
                        eraseOld=FALSE,...
  )
  {
    #browser()
    outval<-compareChoices(assay(x), dimReduce=dimReduce,nVarDims=nVarDims,nPCADims=nPCADims,
                           transFun=transformation(x),...)
    #browser()
    if(class(outval)=="ClusterCells"){
      ##Check if compareChoices already ran previously
      ppIndex<-pipelineClusterIndex(x,print=FALSE)
      if(!is.null(ppIndex)){ #need to change the clusterType values (or erase them) before get new ones
        if(eraseOld){ #remove all of them, not just current
          #browser()
          newX<-removeClusters(x,ppIndex[,"index"]) ###Getting error: Error: evaluation nested too deeply: infinite recursion / options(expressions=)?
        }
        else{
          if(0 %in% ppIndex[,"iteration"]){
            newIteration<-max(ppIndex[,"iteration"])+1
            whCurrent<-ppIndex[ppIndex[,"iteration"]==0,"index"]
            updateCluster<-clusterType(x)
            updateCluster[whCurrent]<-paste(updateCluster[whCurrent],newIteration,sep="_")
            newX<-x
            newX@clusterType<-updateCluster          
          }
        }
        
      }
      else newX<-x
      
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
  definition = function(x, dimReduce="none",nVarDims=NA,nPCADims=NA,
                        transFun=NULL,isCount=FALSE,
                        ...
  )
  {
    outval<-compareChoices(assay(x), 
             dimReduce=dimReduce,nVarDims=nVarDims,nPCADims=nPCADims,
             transFun=transFun,isCount=isCount,...
    )
    if(class(outval)=="ClusterCells"){
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



