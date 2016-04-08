#' @rdname clusterCells
setMethod(
  f = "transform",
  signature = "ClusterCells",
  definition = function(x,nPCADims=NA,nVarDims=NA,dimReduce="none") {
    fun<-transformation(x)
    dat<-assay(x)
    return(.transData(dat,transFun=fun,nPCADims=nPCADims,nVarDims=nVarDims,dimReduce=dimReduce)$x)
  }
)

#function to transform assay data into clustering data (or other normal-like data input)
#if npcs=NA or length of npcs=1, returns matrix; otherwise returns list of pc reduced data.
.transData<-function(x,transFun=NULL,isCount=FALSE,nPCADims,nVarDims,dimReduce)
{
  origX<-x
  #transform data
  if(is.null(transFun)){
    transFun<-if(isCount) function(x){log(x+1)} else function(x){x}
  }
  x<-try(transFun(x),silent=TRUE)
  if(inherits(x, "try-error")) stop(paste("User-supplied `transFun` produces error on the input data matrix:\n",x))
  if(any(is.na(x))) stop("User-supplied `transFun` produces NA values")
  
  ###################
  ###Dim Reduction
  ###################
  ##Check user inputs
  ###################
  if(any(is.na(nPCADims)) & "PCA" %in% dimReduce){
    if(length(nPCADims)==1){
      if(length(dimReduce)==1) dimReduce<-"none" #assume user goofed and meant to do none
      if(length(dimReduce)>1) dimReduce<-dimReduce[-match("PCA",dimReduce)] #assume user goofed and didn't mean to also include 
    } 
    else{
      #add 'none' option to dimReduce and get rid of NA
      dimReduce<-unique(c("none",dimReduce)) 
      nPCADims<-nPCADims[!is.na(nPCADims)]
    }
  }
  if(any(is.na(nVarDims)) & "mostVar" %in% dimReduce){ 
    if(length(nVarDims)==1){
      if(length(dimReduce)==1) dimReduce<-"none" #assume user goofed and meant to do none
      if(length(dimReduce)>1) dimReduce<-dimReduce[-match("mostVar",dimReduce)] #assume user goofed and didn't mean to also include 
      
    } 
    else{#assume user meant to do none as well as dimReduce with other values.
      dimReduce<-unique(c("none",dimReduce)) #add 'none' and remove NA
      nVarDims<-nVarDims[!is.na(nVarDims)]
    }
  }
  
  xPCA<-xVAR<-xNone<-NULL #possible values
  listReturn<-FALSE
  #for each dim reduction method requested
  if("PCA" %in% dimReduce & !all(is.na(nPCADims))){ #do PCA dim reduction
    if(max(nPCADims)>=NROW(x)) stop("the number of PCA dimensions must be strictly less than the number of rows of input data matrix")
    if(min(nPCADims)<1) stop("the number of PCA dimensions must be equal to 1 or greater")
    if(max(nPCADims)>100) warning("the number PCA dimensions to be selected is greater than 100. Are you sure you meant to choose to use PCA dimensionality reduction rather than the top most variable features?")
    prc<-t(stats::prcomp(t(x))$x)
    if(NCOL(prc)!=NCOL(origX)) stop("error in coding of principle components.")
    if(length(nPCADims)==1 & length(dimReduce)==1){ #just return single matrix
      x<-prc[1:nPCADims,]
      
    }
    else{
      xPCA<-lapply(nPCADims,function(nn){prc[1:nn,]})
      names(xPCA)<-paste("nPCAFeatures=",nPCADims,sep="")
      listReturn<-TRUE
    }
  }
  if("mostVar" %in% dimReduce & all(!is.na(nVarDims))){ #do PCA dim reduction
    if(max(nVarDims)>=NROW(x)) stop("the number of most variable features must be strictly less than the number of rows of input data matrix")
    if(min(nVarDims)<1) stop("the number of most variable features must be equal to 1 or greater")
    if(min(nVarDims)<50 & NROW(x)>1000) warning("the number of most variable features to be selected is less than 50. Are you sure you meant to choose to use the top most variable features rather than PCA dimensionality reduction?")
    varX<-apply(x,1,mad)
    ord<-order(varX,decreasing=TRUE)
    xVarOrdered<-x[ord,]
    if(NCOL(xVarOrdered)!=NCOL(origX)) stop("error in coding of principle components.")
    if(length(nVarDims)==1 & length(dimReduce)==1){ #just return single matrix
      x<-xVarOrdered[1:nVarDims,]
    }
    else{ #otherwise make it a list
      xVAR<-lapply(nVarDims,function(nn){xVarOrdered[1:nn,]})
      names(xVAR)<-paste("nVarFeatures=",nVarDims,sep="")
      listReturn<-TRUE
    }
  }
  if("none" %in% dimReduce & length(dimReduce)>1){
    xNone<-list("noDimReduce"=x)
    listReturn<-TRUE
  }
  #browser()
  
  if(listReturn) x<-c(xNone,xVAR,xPCA)
  return(list(x=x,transFun=transFun))
}