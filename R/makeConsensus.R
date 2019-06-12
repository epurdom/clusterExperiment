#'Find sets of samples that stay together across clusterings
#'
#'Find sets of samples that stay together across clusterings in order to define 
#'a new clustering vector.
#'
#' @aliases makeConsensus
#'  
#' @param x a matrix or \code{\link{ClusterExperiment}} object.

#' @param clusterFunction the clustering to use (passed to 
#'  \code{\link{mainClustering}}); currently must be of type '01'.
#' @param minSize minimum size required for a set of samples to be considered in 
#'  a cluster because of shared clustering, passed to
#'  \code{\link{mainClustering}}
#' @param proportion The proportion of times that two sets of samples should be 
#'  together in order to be grouped into a cluster (if <1, passed to
#'  mainClustering via alpha = 1 - proportion)
#' @param propUnassigned samples with greater than this proportion of assignments
#'  equal to '-1' are assigned a '-1' cluster value as a last step (only if
#'  proportion < 1)
#' @param ... arguments to be passed on to the method for signature 
#'  \code{matrix,missing}.
#' @inheritParams clusterMany
#' @inheritParams getClusterIndex
#' @details This function was previously called \code{combineMany} (versions <=
#'   2.0.0). \code{combineMany} is still available, but is considered defunct
#'   and users should update their code accordingly.
#' @details The function tries to find a consensus cluster across many different 
#'  clusterings of the same samples. It does so by creating a \code{nSamples} x 
#'  \code{nSamples} matrix of the percentage of co-occurance of each sample and 
#'  then calling mainClustering to cluster the co-occurance matrix. The function
#'  assumes that '-1' labels indicate clusters that are not assigned to a 
#'  cluster. Co-occurance with the unassigned cluster is treated differently 
#'  than other clusters. The percent co-occurance is taken only with respect to 
#'  those clusterings where both samples were assigned. Then samples with more 
#'  than \code{propUnassigned} values that are '-1' across all of the 
#'  clusterings are assigned a '-1' regardless of their cluster assignment.
#'@details The method calls \code{\link{mainClustering}} on the proportion
#'  matrix with \code{clusterFunction} as the 01 clustering algorithm,
#'  \code{alpha=1-proportion}, \code{minSize=minSize}, and
#'  \code{evalClusterMethod=c("average")}. See help of 
#'  \code{\link{mainClustering}} for more details.
#'@return If x is a matrix, a list with values \itemize{ 
#'  \item{\code{clustering}}{ vector of cluster assignments, with "-1" implying 
#'  unassigned}
#'  
#'  \item{\code{percentageShared}}{ a nSample x nSample matrix of the percent 
#'  co-occurance across clusters used to find the final clusters. Percentage is 
#'  out of those not '-1'} \item{\code{noUnassignedCorrection}{ a vector of 
#'  cluster assignments before samples were converted to '-1' because had 
#'  >\code{propUnassigned} '-1' values (i.e. the direct output of the 
#'  \code{mainClustering} output.)}} }
#'  
#' @return If x is a \code{\link{ClusterExperiment}}, a
#'  \code{\link{ClusterExperiment}} object, with an added clustering of
#'  clusterTypes equal to \code{makeConsensus} and the \code{percentageShared}
#'  matrix stored in the \code{coClustering} slot.
#'
#' @examples
#' data(simData)
#'
#' cl <- clusterMany(simData,nReducedDims=c(5,10,50),  reduceMethod="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(FALSE), removeSil=TRUE,
#' subsample=FALSE)
#'
#' #make names shorter for plotting
#' clMat <- clusterMatrix(cl)
#' colnames(clMat) <- gsub("TRUE", "T", colnames(clMat))
#' colnames(clMat) <- gsub("FALSE", "F", colnames(clMat))
#' colnames(clMat) <- gsub("k=NA,", "", colnames(clMat))
#'
#' #require 100% agreement -- very strict
#' clCommon100 <- makeConsensus(clMat, proportion=1, minSize=10)
#'
#' #require 70% agreement based on clustering of overlap
#' clCommon70 <- makeConsensus(clMat, proportion=0.7, minSize=10)
#'
#' oldpar <- par()
#' par(mar=c(1.1, 12.1, 1.1, 1.1))
#' plotClusters(cbind("70%Similarity"=clCommon70, clMat,
#' "100%Similarity"=clCommon100), axisLine=-2)
#'
#' #method for ClusterExperiment object
#' clCommon <- makeConsensus(cl, whichClusters="workflow", proportion=0.7,
#' minSize=10)
#' plotClusters(clCommon)
#' par(oldpar)
#'
#' @rdname makeConsensus
#' @export
setMethod(
    f = "makeConsensus",
    signature = signature(x = "matrix"),
    definition = function(x, proportion,
                          clusterFunction="hierarchical01",
                          minSize=5, propUnassigned=.5, 
                          whenUnassign=c("before","after","never"),
                          clusterArgs=NULL) {
        whenUnassign<-match.arg(whenUnassign)
        if(proportion >1 || proportion <0) stop("Invalid value for the 'proportion' parameter")
        if(propUnassigned >1 || propUnassigned <0) stop("Invalid value for the 'propUnassigned' parameter")
        #FIXME: unnecessary duplication
        clusterMat <- x
        ##Now define as unassigned any samples with >= propUnassigned '-1' values in clusterMat
        whUnassigned <- which(apply(clusterMat, 2, function(x){
            sum(x== -1)/length(x)>propUnassigned}))
        N<-nrow(clusterMat)
        if(length(whUnassigned)>0 && whUnassigned=="before") clusterMat<-clusterMat[-whUnassigned,]
        if(proportion == 1) {
            #have to repeat from mainClustering because didn't
            if(!is.numeric(minSize) || minSize<0) 
                stop("Invalid value for the 'minSize' parameter in determining the minimum number of samples required in a cluster.")
            else minSize<-round(minSize) #incase not integer.
            singleValueClusters <- apply(clusterMat, 1, paste, collapse=";")
            allUnass <- paste(rep("-1", length=ncol(clusterMat)), collapse=";")
            uniqueSingleValueClusters <- unique(singleValueClusters)
            tab <-	table(singleValueClusters)
            tab <- tab[tab >= minSize]
            tab <- tab[names(tab) != allUnass]
            cl <- match(singleValueClusters, names(tab))
            cl[is.na(cl)] <- -1
        } else{
            
            if(is.character(clusterFunction)) 
                typeAlg <- algorithmType(clusterFunction)
            else{
                if(is(clusterFunction,"ClusterFunction")) 
                    typeAlg<-algorithmType(clusterFunction) 
                else 
                    stop("clusterFunction must be either built in clusterFunction name or a ClusterFunction object")
            }
            if(typeAlg!="01") {
                stop("makeConsensus is only implemented for '01' type clustering functions (see ?ClusterFunction)")
            }
            #overwrite if given
            clusterArgs[["alpha"]]<-1-proportion
            if(!"evalClusterMethod" %in% names(clusterArgs) &&
                clusterFunction=="hierarchical01"){
                clusterArgs<-c(clusterArgs,
                    list(evalClusterMethod=c("average")))
            }
            cl <- mainClustering(inputMatrix=t(clusterMat),
                                 inputType="cat",
                                 clusterFunction=clusterFunction,
                                 minSize=minSize, format="vector",
                                 clusterArgs=clusterArgs)
            
            if(is.character(cl)) {
                stop("coding error -- mainClustering should return numeric vector")
            }
        }
        if(length(whUnassigned)>0){
            if(whUnassigned=="after") cl[whUnassigned]<- -1
            if(whUnassigned=="before"){
                #put back in the -1 values
                temp<-rep(-1,length=N)
                if(length(whUnassigned)!=N) temp[-whUnassigned]<-cl
                cl<-temp
                rm(temp)
            }                    
        }
        return(cl)
    }
)

.unassignLargeMissing<-function(clusterMat,propUnassigned){
    
}

#' @rdname makeConsensus
#' @export
#' @param clusterLabel a string used to describe the type of clustering. By
#'   default it is equal to "makeConsensus", to indicate that this clustering is
#'   the result of a call to makeConsensus. However, a more informative label can
#'   be set (see vignette).
setMethod(
    f = "makeConsensus",
    signature = signature(x = "ClusterExperiment"),
    definition = function(x, whichClusters, eraseOld=FALSE,clusterLabel="makeConsensus",...){
        if(missing(whichClusters)){
            whichClusters <- getClusterIndex(x, whichClusters="clusterMany", noMatch="silentlyRemove")
            if(length(whichClusters)>0){
                .mynote("no clusters specified to combine, using results from clusterMany")
            }
            else{
                stop("no clusters specified to combine, please specify.")
            }
        }
        else{
            whichClusters <-getClusterIndex(x,whichClusters=whichClusters,noMatch="throwError")
        }
        outlist <- makeConsensus(clusterMatrix(x)[, whichClusters, drop=FALSE], ...)
        newObj <- ClusterExperiment(x, outlist,
                                    transformation=transformation(x),
                                    clusterTypes="makeConsensus",checkTransformAndAssay=FALSE)
        #add "c" to name of cluster
        newObj<-.addPrefixToClusterNames(newObj,prefix="c",whCluster=1)
        clusterLabels(newObj) <- clusterLabel
        
        # ## FIXME: the function no longer returns this object.
        #  need to the index of clusterings used for the makeConsensus?
        # Or save it only if plot it (i.e. calculate it)?
        # if(!is.null(outlist$percentageShared)) {
        #   coClustering(newObj) <- Matrix::Matrix(outlist$percentageShared,sparse=TRUE)
        # }
        coClustering(newObj) <- clusterMatrix(x)[, whichClusters, drop=FALSE]
        ##Check if pipeline already ran previously and if so increase
        x<-.updateCurrentWorkflow(x,eraseOld,newTypeToAdd="makeConsensus",newLabelToAdd=clusterLabel)
        
        if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x) #make decisions about what to keep.
        else retval<-.addBackSEInfo(newObj=newObj,oldObj=x)
        return(retval)
    }
)




