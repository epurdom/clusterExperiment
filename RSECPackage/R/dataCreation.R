#' @title Documentation of rsecFluidigm object
#' @description  Documentation of the creation of rsecFluidigm, result of RSEC
#'   run on fluidigm data for vignette
#'
#' @name rsecFluidigm
#' @docType data
#' @author Elizabeth Purdom \email{epurdom@@stat.berkeley.edu}
#' @format \code{rsecFluidigm} is a \code{ClusterExperiment} object, the result
#'   of running \code{\link{RSEC}} on fluidigm data described in vignette and
#'   available in the \code{scRNAseq} package.
#' @details The functions \code{makeRsecFluidigmObject} and 
#'   \code{checkRsecFluidigmObject} are helper functions whose sole purpose is 
#'   to create \code{rsecFluidigm} and check that the results are the same as 
#'   expected. \code{makeRsecFluidigmObject} also serves as documentation of the
#'   specific RSEC call that was made to create the \code{rsecFluidigm} object, 
#'   as well as filtering and normalization of the fluidigm data.
#'   The purpose of making them functions is internal, to help more easily
#'   mantain and check if changes to the package have affected the results.
#' @seealso \code{\link[scRNAseq]{fluidigm}}. 
#' @keywords data
#' @examples
#' # see code used create rsecFluidigm 
#' # (print out the function)
#' makeRsecFluidigmObject
#' #code actualy run to create rsecFluidigm:
#' \dontrun{
#' library(clusterExperiment)
#' data(fluidigmData)
#' data(fluidigmColData)
#' se<-SummarizedExperiment(assays=fluidigmData, colData=fluidigmColData)
#' RNGversion("3.5.0")
#' rsecFluidigm<-makeRsecFluidigmObject(se)
#' # Internal function for checking got correct results...
#' clusterExperiment:::checkRsecFluidigmObject(rsecFluidigm)
#' usethis::use_data(rsecFluidigm,overwrite=FALSE)
#' }
#' @aliases makeRsecFluidigmObject
#' @param object object given to functions
#' @export
#' @importFrom S4Vectors metadata
makeRsecFluidigmObject<-function(object){
    pass_filter <- apply(assay(object), 1, 
        function(x) length(x[x >= 10]) >= 10)
    object <- object[pass_filter,]
    fq <- round(limma::normalizeQuantiles(assay(object)))
    assays(object) <- c(SimpleList(normalized_counts=fq),assays(object))
    wh<-which(colnames(colData(object)) %in% c("Cluster1","Cluster2"))
    colnames(colData(object))[wh]<-c("Published1","Published2")
    ncores<-1
    rsecFluidigm<-RSEC(object,
                      isCount = TRUE,
                      k0s = 4:15,
                      alphas=c(0.1, 0.2, 0.3),
                      betas = 0.9,
                      reduceMethod="PCA",
                      nReducedDims=10,
                      minSizes=1,
                      subsample=TRUE,
                      sequential=TRUE,
                      clusterFunction="hierarchical01",
                      consensusMinSize=3,
                      consensusProportion=0.7,
                      dendroReduce= "mad",
                      dendroNDims=1000,
                      mergeMethod="adjP",
                      mergeDEMethod="limma",
                      mergeCutoff=0.01,
                      ncores=ncores,
                      makeMissingDiss=TRUE,
                      mainClusterArgs=list(
                          clusterArgs=list(removeDup=FALSE)),
                      #seqArgs=list(top.can=5),
                      subsampleArgs=list(clusterFunction="kmeans",
                          classifyMethod="All"),
                      consensusArgs=list(clusterFunction="hierarchical01",
                            whenUnassign="before",
                            clusterArgs=list(
                                evalClusterMethod=c("average"),
                                removeDup=FALSE)),
                      random.seed=176201
    )
    metadata(rsecFluidigm)$packageVersion <- packageVersion("clusterExperiment")
    return(rsecFluidigm)
}

checkRsecFluidigmObject<-function(object){
    ## Simple Tests that haven't changed the clustering algorithms such that get different results.
    ## Don't simply do all.equal with old one because might of changed something minor not related to the actual algorithms

    ## Results for feature/knn 07/25/2019 -- 2.5.4.9005
    nMakeConsensus<-7
    nMerge<-5
    contrasts<-c('(X4+X1+X5)/3-(X6+X7+X2+X3)/4','X4-(X1+X5)/2',
        'X6-(X7+X2+X3)/3','X7-(X2+X3)/2','X1-X5','X2-X3')                  
    adjPValues<-c(0.06167775,0.01117556,0.01697553,0.00042439,0.01004385,0.00155609)
       
    ## Test same
    checkValues<-.getCheckValues(object)
    if(nMakeConsensus==nMerge) warning("having same number of clusters before and after merge won't be great for the vignette!")
    if(checkValues$nMakeConsensus != nMakeConsensus)
      stop("rsecFluidigm has changed -- makeConsensus")
    if(!is.logical(all.equal(contrasts,checkValues$contrasts)))
        stop("rsecFluidigm has changed -- different dendrogram")
    if(!is.logical(all.equal(adjPValues,checkValues$adjPValues)))
      stop("rsecFluidigm has changed -- different merge percentages")
    if(checkValues$nMerge != nMerge)
      stop("rsecFluidigm has changed -- different # of mergeClusters")
    
}
.getCheckValues<-function(object, printout=FALSE){
    x<-unique(clusterMatrix(object)[,"makeConsensus"])
    y<-unique(clusterMatrix(object)[,"mergeClusters"])
    out<-list(
        nMakeConsensus=length(x[x>0]),
        nMerge=length(y[y>0]),
        contrasts=object@merge_nodeProp[,"Contrast"],
        adjPValues=round(object@merge_nodeProp[,"adjP"],8)
        )
    if(printout){
        out$contrasts<-paste("c('",paste(out$contrasts,collapse="','"),"')",sep="")
        out$adjPValues<-paste("c(",paste(out$adjPValues,collapse=","),")",sep="")
    }
    return(out)
}





#' Subset of fluidigm data
#'
#' @name fluidigmData
#' @aliases fluidigmColData
#' @docType data
#' @author Elizabeth Purdom \email{epurdom@@stat.berkeley.edu}
#' @format subset of fluidigm data used in vignette
#' package.
#' @seealso \code{\link[scRNAseq]{fluidigm}}
#' @keywords data
#' @details \code{fluidigmData} and \code{fluidigmColData} are portions of the \code{fluidigm} data distributed in the package \code{scRNAseq} package. We have subsetted to only the cells sequenced under high depth, and limited our selves to only two of the four gene estimates provided by \code{scRNAseq} ("tophat_counts" and "rsem_tpm").
#' @examples
#' #code used to create objects:
#' \dontrun{
#' library(scRNAseq)
#' if(packageVersion("scRNAseq")>="1.11.0") fluidigm <- ReprocessedFluidigmData() else data(fluidigm)
#' fluidSubset<- fluidigm[,colData(fluidigm)[,"Coverage_Type"]=="High"]
#' fluidigmData<-assays(fluidSubset)[c("tophat_counts","rsem_tpm")]
#' fluidigmColData<-as.data.frame(colData(fluidSubset))
#' usethis::use_data(fluidigmData, fluidigmColData, overwrite=FALSE)
#' }
NULL
