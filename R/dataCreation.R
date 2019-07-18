#' Simulated data for running examples
#'
#' @name simData
#' @aliases simCount trueCluster
#' @docType data
#' @author Elizabeth Purdom \email{epurdom@@stat.berkeley.edu}
#' @format Three objects are loaded, two data frame(s) of simulated data each
#'   with 300 samples/columns and 153 variables/rows, and a vector of length 300
#'   with the true cluster assignments.
#' @details \code{simData} is simulated normal data of 300 observations with 51
#'   relevant variables and the rest of the variables being noise, with
#'   observations being in one of 3 groups. \code{simCount} is simulated count
#'   data of the same dimensions. \code{trueCluster} gives the true cluster
#'   identifications of the samples. The true clusters are each of size 100 and
#'   are in order in the columns of the data.frames.
#' @keywords data
#' @examples
#' #code used to create data:
#' \dontrun{
#' nvar<-51 #multiple of 3
#' n<-100
#' x<-cbind(matrix(rnorm(n*nvar,mean=5),nrow=nvar),
#'  matrix(rnorm(n*nvar,mean=-5),nrow=nvar),
#'           matrix(rnorm(n*nvar,mean=0),nrow=nvar))
#' #make some of them flipped effects (better for testing if both sig under/over
#' #expressed variables)
#' geneGroup<-sample(rep(1:3,each=floor(nvar/3)))
#' gpIndex<-list(1:n,(n+1):(n*2),(2*n+1):(n*3))
#' x[geneGroup==1,]<-x[geneGroup==1,unlist(gpIndex[c(3,1,2)])]
#' x[geneGroup==2,]<-x[geneGroup==2,unlist(gpIndex[c(2,3,1)])]
#'
#' #add in differences in variable means
#' smp<-sample(1:nrow(x),10)
#' x[smp,]<-x[smp,]+10
#'
#' #make different signal y
#' y<-cbind(matrix(rnorm(n*nvar,mean=1),nrow=nvar),
#'          matrix(rnorm(n*nvar,mean=-1),nrow=nvar),
#'          matrix(rnorm(n*nvar,mean=0),nrow=nvar))
#' y<-y[,sample(1:ncol(y))]+ matrix(rnorm(3*n*nvar,sd=3),nrow=nvar)
#'
#' #add together the two signals
#' simData<-x+y
#'
#' #add pure noise variables
#' simData<-rbind(simData,matrix(rnorm(3*n*nvar,mean=10),nrow=nvar),
#'                matrix(rnorm(3*n*nvar,mean=5),nrow=nvar))
#' #make count data
#' countMean<-exp(simData/2)
#' simCount<-matrix(rpois(n=length(as.vector(countMean)), lambda
#' =as.vector(countMean)+.1),nrow=nrow(countMean),ncol=ncol(countMean))
#' #labels for the truth
#' trueCluster<-rep(c(1:3),each=n)
#' save(list=c("simCount","simData","trueCluster"),file="data/simData.rda")
#' }
NULL


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
#' se<-SummarizedExperiment(assays=fluidigmobjectolData=fluidigmColData)
#' rsecFluidigm<-makeRsecFluidigmObject(se)
#' checkRsecFluidigmObject(rsecFluidigm)
#' usethis::use_data(rsecFluidigm,overwrite=FALSE)
#' }
#' @aliases makeRsecFluidigmObject
#' @param object object given to functions
#' @export
makeRsecFluidigmObject<-function(object){
    pass_filter <- apply(assay(object), 1, function(x) length(x[x >= 10]) >= 10)
    object <- object[pass_filter,]
    fq <- round(limma::normalizeQuantiles(assay(object)))
    assays(object) <- c(SimpleList(normalized_counts=fq),assays(object))
    wh<-which(colnames(colData(object)) %in% c("Cluster1","Cluster2"))
    colnames(colData(object))[wh]<-c("Published1","Published2")
    ncores<-1
    system.time(
      rsecFluidigm<-RSEC(object, 
                         isCount = TRUE, 
                         k0s = 4:15, 
                         alphas=c(0.1, 0.2, 0.3), 
                         betas = 0.9,
                         reduceMethod="PCA", 
                         nReducedDims=10,
                         minSizes=1, 
                         clusterFunction="hierarchical01",
                         consensusMinSize=3,
                         consensusProportion=0.7,
                         dendroReduce= "mad",
                         dendroNDims=1000,
                         mergeMethod="adjP",
    	                    mergeDEMethod="limma",
                         mergeCutoff=0.01,
                         ncores=ncores, 
                         random.seed=176201)
    )
    SummarizedExperiment::metadata(rsecFluidigm)$packageVersion <- packageVersion("clusterExperiment")
    return(rsecFluidigm)
}
#' @rdname makeRsecFluidigmObject
#' @export
checkRsecFluidigmObject<-function(object){
    ## Simple Tests that haven't changed the clustering algorithms such that get different results.
    ## Don't simply do all.equal with old one because might of changed something minor not related to the actual algorithms

    ## Expected Results
    
    ## Previous expected results -- can't recreate!
    # nMakeConsensus<-8
    # nMerge<-6
    # contrasts<-c('(X6+X1+X4)/3-(X8+X2+X5+X3+X7)/5','X6-(X1+X4)/2','X8-(X2+X5+X3+X7)/4','(X2+X5)/2-(X3+X7)/2','X2-X5','X3-X7','X1-X4')
    # adjPValues<-c(0.049794879, 0.007356062, 0.008204838,
    #               0.013156033, 0.009336540, 0.007497524, 0.033526666)

    ## Results as of 07/12/2019 -- 2.5.4.9002
    nMakeConsensus<-10
    nMerge<-6
    contrasts<-c('(X10+X2+X3+X6+X8)/5-(X4+X5+X9+X1+X7)/5','(X10+X2+X3)/3-(X6+X8)/2','(X4+X5)/2-(X9+X1+X7)/3','X4-X5','X9-(X1+X7)/2','X1-X7','X10-(X2+X3)/2','X6-X8','X2-X3')
    adjPValues<-c(0.08374593,0.03762908,0.01372188,0.0072146,0.00778045,0.01018532,0.00056585,0.00650729,0.00183902)

        
    ## Test same
    checkValues<-.getCheckValues(object)
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


