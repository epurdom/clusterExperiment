#' Program for sequentially clustering, removing cluster, and starting again.
#'
#' Given a data matrix, this function will call clustering
#' routines, and sequentially remove best clusters, and iterate to find
#' clusters.
#'
#' @param k0 the value of K at the first iteration of sequential algorithm, see
#'   details below or vignette.
#' @param subsample logical as to whether to subsample via 
#'   \code{\link{subsampleClustering}} to get the distance matrix at each 
#'   iteration; otherwise the distance matrix is set by arguments to
#'   \code{\link{mainClustering}}.
#' @param beta value between 0 and 1 to decide how stable clustership membership
#'   has to be before 'finding' and removing the cluster.
#' @param top.can only the top.can clusters from \code{\link{mainClustering}} (ranked
#'   by 'orderBy' argument given to \code{\link{mainClustering}}) will be compared
#'   pairwise for stability. Can be either an integer value, identifying the absolute number of clusters, or a value between 0 and 1, meaning to keep all clusters with at least this proportion of the remaining samples in the cluster. Making this either a very big integer or equal to 0 will effectively remove this
#'   parameter and all pairwise comparisons of all clusters found will be
#'   considered; this might result in smaller clusters being found. If top.can is between 0 and 1, then there is still a hard threshold of at least 5 samples in a cluster to be considered as a cluster. 
#' @param remain.n when only this number of samples are left (i.e. not yet
#'   clustered) then algorithm will stop.
#' @param k.min each iteration of sequential detection of clustering will
#'   decrease the beginning K of subsampling, but not lower than k.min. 
#' @param k.max algorithm will stop if K in iteration is increased beyond this
#'   point.
#' @param verbose whether the algorithm should print out information as to its
#'   progress.
#' @param subsampleArgs list of arguments to be passed to
#'   \code{\link{subsampleClustering}}.
#' @param mainClusterArgs list of arguments to be passed to
#'   \code{\link{mainClustering}}).
#' @inheritParams clusterSingle
#' @details \code{seqCluster} is not meant to be called by the user. It is only
#'   an exported function so as to be able to clearly document the arguments for
#'   \code{seqCluster} which can be passed via the argument \code{seqArgs} in
#'   functions like \code{\link{clusterSingle}} and \code{\link{clusterMany}}.
#' @details This code is adapted from the sequential protion of the code of the
#'   tightClust package of Tseng and Wong. At each iteration of the algorithm it
#'   finds a set of samples that constitute a homogeneous cluster and remove
#'   them, and iterate again to find the next set of samples that form a
#'   cluster.
#' @details In each iteration, to determine the next set of homogeneous set of
#'   samples, the algorithm will iteratively cluster the current set of samples
#'   for a series of increasing values of the parameter $K$, starting at a value
#'   \code{kinit} and increasing by 1 at each iteration, until a sufficiently
#'   homogeneous set of clusters is found. For the first set of homogeneous
#'   samples, \code{kinit} is set to the argument $k0$, and for iteration,
#'   \code{kinit} is increased internally.
#' @details Depending on the value of \code{subsample} how the value of $K$ is
#'   used differs. If \code{subsample=TRUE}, $K$ is the \code{k} sent to the
#'   cluster function \code{clusterFunction} sent to 
#'   \code{\link{subsampleClustering}} via \code{subsampleArgs}; then
#'   \code{\link{mainClustering}} is run on the result of the co-occurance matrix from
#'   \code{\link{subsampleClustering}} with the \code{ClusterFunction} object
#'   defined in the argument \code{clusterFunction} set via \code{mainClusterArgs}.
#'   The number of clusters actually resulting from this run of
#'   \code{\link{mainClustering}} may not be equal to the $K$ sent to  the clustering
#'   done in \code{\link{subsampleClustering}}. If \code{subsample=FALSE},
#'   \code{\link{mainClustering}} is called directly on the data to determine the
#'   clusters and $K$ set by \code{seqCluster} for this iteration determines the
#'   parameter of the clustering done by \code{\link{mainClustering}}. Specifically,
#'   the argument \code{clusterFunction} defines the clustering of the
#'   \code{\link{mainClustering}} step and \code{k} is sent to that
#'   \code{ClusterFunction} object. This means that if \code{subsample=FALSE},
#'   the \code{clusterFunction} must be of \code{algorithmType} "K".
#' @details In either setting of \code{subsample}, the resulting clusters from
#'   \code{\link{mainClustering}} for a particular $K$ will be compared to clusters
#'   found in the previous iteration of $K-1$. For computational (and other?)
#'   convenience, only the first \code{top.can} clusters of each iteration will
#'   be compared to the first \code{top.can} clusters of previous iteration for
#'   similarity (where \code{top.can} currently refers to ordering by size, so
#'   first \code{top.can} largest clusters.
#' @details If there is no cluster of the first \code{top.can} in the current
#'   iteration $K$ that has overlap similarity > \code{beta} to any in the
#'   previous iteration, then the algorithm will move to the next iteration,
#'   increasing to $K+1$.
#'   
#' @details If, however, of these clusters there is a cluster in the current
#'   iteration $K$ that has overlap similarity > beta to a cluster in the
#'   previous iteration $K-1$, then the cluster with the largest such similarity
#'   will be identified as a homogenous set of samples and the samples in it
#'   will be removed and designated as such. The algorithm will then start again
#'   to determine the next set of homogenous samples, but without these samples.
#'   Furthermore, in this case (i.e. a cluster was found and removed), the value
#'   of \code{kinit} will be be reset to \code{kinit-1}; i.e. the range of
#'   increasing $K$ that will be iterated over to find a set of homogenous
#'   samples will start off one value less than was the case for the previous
#'   set of homogeneous samples. If \code{kinit-1}<\code{k.min}, then
#'   \code{kinit} will be set to \code{k.min}.
#'   
#'   
#' @details If there are less than \code{remain.n} samples left after finding a
#'   cluster and removing its samples, the algorithm will stop, as subsampling
#'   is deamed to no longer be appropriate. If the K has to be increased to
#'   beyond \code{k.max} without finding any pair of clusters with overlap >
#'   beta, then the algorithm will stop. Any samples not found as part of a
#'   homogenous set of clusters at that point will be classified as unclustered
#'   (given a value of -1)
#' @details Certain combinations of inputs to \code{mainClusterArgs} and
#'   \code{subsampleArgs} are not allowed. See \code{\link{clusterSingle}} for
#'   these explanations.
#' @return A list with values
#' \itemize{
#'
#' \item{\code{clustering}}{ a vector of length equal to nrows(x) giving the
#' integer-valued cluster ids for each sample. The integer values are assigned
#' in the order that the clusters were found. "-1" indicates the sample was not
#' clustered.}
#'
#' \item{\code{clusterInfo}}{ if clusters were successfully found, a matrix of
#' information regarding the algorithm behavior for each cluster (the starting
#' and stopping K for each cluster, and the number of iterations for each
#' cluster).}
#'
#' \item{\code{whyStop}}{ a character string explaining what triggered the
#' algorithm to stop.}
#' }
#' @references Tseng and Wong (2005), "Tight Clustering: A Resampling-Based
#'   Approach for Identifying Stable and Tight Patterns in Data", Biometrics,
#'   61:10-16.
#' 
#' @seealso tight.clust,
#'   \code{\link{clusterSingle}},\code{\link{mainClustering}},\code{\link{subsampleClustering}}
#'   
#' @examples
#' \dontrun{
#' data(simData)
#'
#' set.seed(12908)
#' clustSeqHier <- seqCluster(simData, inputType="X", k0=5, subsample=TRUE,
#'    beta=0.8, subsampleArgs=list(resamp.n=100,
#'    samp.p=0.7, clusterFunction="kmeans", clusterArgs=list(nstart=10)),
#'    mainClusterArgs=list(minSize=5,clusterFunction="hierarchical01",
#'    clusterArgs=list(alpha=0.1)))
#' }
#' @export
#' @rdname seqCluster
#' @export
seqCluster<-function(inputMatrix, inputType, k0,  
                     subsample=TRUE, beta, 
                     top.can = 0.01, remain.n = 30, k.min = 3, 
                     k.max=k0+10,verbose=TRUE, 
                     subsampleArgs=NULL,mainClusterArgs=NULL,
                     warnings=FALSE)
{
    ########
    ####Checks
    ########
    if(missing(inputType)) 
        stop("Internal error: inputType was not passed to sequential clustering step")
    argsList<-list(k0=k0, beta=beta, top.can = top.can, 
        remain.n = remain.n, k.min = k.min, k.max=k.max)
    checkOut<-.checkArgs(inputType, main=TRUE, 
        subsample=subsample,sequential=TRUE,
		mainClusterArgs=mainClusterArgs,
		subsampleArgs=subsampleArgs, seqArgs=argsList,warn=warnings)
    if(is.character(checkOut)) stop(checkOut)
    mainClusterArgs<-checkOut$mainClusterArgs
    subsampleArgs<-checkOut$subsampleArgs
    seqArgs<-checkOut$seqArgs
    k0=seqArgs[["k0"]]
    beta=seqArgs[["beta"]]
    top.can = seqArgs[["top.can"]]
    isIntegerTop<-is.whole(top.can) & top.can >0
    remain.n = seqArgs[["remain.n"]]
    k.min = seqArgs[["k.min"]]
    k.max=seqArgs[["k.max"]]
    
    if(subsample) input<-subsampleArgs[["inputType"]]
    else input<-mainClusterArgs[["inputType"]]
    #Note, for the purpose of this
    N <- dim(inputMatrix)[2]
    colnames(inputMatrix)<-as.character(seq_len(N))
    if(input=="diss") rownames(inputMatrix)<-colnames(inputMatrix)
    
    #########
    # Function that updates k and runs clustering
    # (Called at each iteration)
    # Note that depends on .checkArgs having removed the 'k' option 
    # from the user inputs
    #########
    updateClustering<-function(newk){
        if(verbose) cat(paste("k =", newk,"\n"))
        # set k to the new k in main cluster step
        # note is subsample=TRUE, and alg. of main is type K, 
        # means k is always equal at subsample and main step. 
        tempMainArgs<-mainClusterArgs
        if(algorithmType(tempMainArgs[["clusterFunction"]])=="K")
            tempMainArgs[["clusterArgs"]]<-c(list(k=newk),
			mainClusterArgs[["clusterArgs"]])       
        else tempMainArgs<-mainClusterArgs
        tempSubArgs<-subsampleArgs
        # set k to new k in subsampling step.
        if(subsample){
            tempSubArgs[["clusterArgs"]]<-c(list(k=newk), subsampleArgs[["clusterArgs"]]) 
        }
        res <- .clusterWrapper(inputMatrix=inputMatrix, 
                subsample=subsample,  subsampleArgs=tempSubArgs,
                mainClusterArgs=tempMainArgs)$results
        return(res)
    }
    trimResults<-function(res){
        if(length(res)>0){
            if(isIntegerTop){
                res <- res[seq_len(min(top.can,length(res)))]
            } 
            else{
                currN<-ncol(inputMatrix)
                requiredN<-floor(currN*top.can)
                NPerCluster<-sapply(res,length)
                if(isTRUE(all.equal(top.can,0))){
                    requiredN<-0
                }
                else{
                    requiredN<-max(requiredN,5)
                }
                whLarger<-which(NPerCluster>=requiredN)
                # cat("Npercluster:\n")
                # print(NPerCluster)
                if(length(whLarger)>0){
                    res<-res[whLarger]
                    #cat(sprintf("Npercluster that larger than %s:\n",max(requiredN,5)))
                    #print(NPerCluster[whLarger])
                }
                else res<-res[seq_len(min(10,length(res)))]
            }
        } 
        return(res)
    }
    #########
    ## Iteration code
    #########
    
    ###---------
    ###---------
    ###The following is legacy of tight.clust. They originally had programmed ability to look across more than 2 at each step to determing the stability of a cluster. This was not what they described in paper, and function is hard-coded at 2, but I have left code here in case we ever wanted to reconsider this issue.
    ###---------
    seq.num<-2
    kReturn<-"last" # when look at stability, return stable as first or last? For seq.num=2, not really matter, take last like paper
    kReturn<-match.arg(kReturn,c("last","first"))
    betaNum<-"all"
    betaNum<-match.arg(betaNum,c("all","last","first"))
    #This makes all combinations of 1:top.can, seq.num times (could be simplified if seq.num=2):
    #a ncombinations x seq.num matrix -- each row gives a combination of clusters to compare stability
    if( isIntegerTop) index.m <- as.matrix(expand.grid(lapply(seq_len(seq.num), function(x) seq_len(top.can))))
    else index.m<-NULL
    whReturn<-switch(kReturn,"last"=seq.num,"first"=1) #way to index which one gets returned.
    ###---------
    ###---------
    
    #---------
    #iterative setup
    #-------
    remain <- N #keep track of how many samples not yet clustered (stop when less than remain.n)
    nfound <- 0 #keep track of how many clusters found/removed so far
    found <- TRUE #has a cluster been found/removed in last iteration
    k.start <- k0 #the starting k for the next cluster
    k <- k0
    ### Blank values that will be filled in.
    candidates <- list() #list of length seq.num of possible clusters found for each k to be compared
    tclust <- list() #list of final cluster identifications (indices of rows of x)
    kstart<-c() #the starting k for the cluster
    kend<-c() #the ending k for the cluster
    whyStop<-NULL
    #---------
    #iterative loop
    #-------
    while (remain >= remain.n && (found || k <= k.max)) {
        if (found) { 
            # i.e. previous iteration found cluster; 
            # need to start finding new cluster
            if(verbose) cat(paste("Looking for cluster", nfound + 1, "...\n"))
            k <- k.start
            currentStart<-k.start #will add this to kstart if successful in finding cluster
            # find clusters for K,K+1
            # Note, candidates should be empty, so filling first two slots
            for (ii in seq_len(seq.num)) {
                newk<-k + ii - 1
                res<-updateClustering(newk)
                candidates[[ii]]<-trimResults(res)
            }
        }
        else { 
            # previous iteration clustering wasn't good enough
            # need to go increase to K+2,K+3, etc.
            candidates <- candidates[-1] #remove old k
            if(length(candidates)!=1 ) stop("Internal coding error in seqCluster -- more candidates than considering")
            newk<-k + seq.num - 1
            #add new k (because always list o)
            res<-updateClustering(newk)
            candidates[[seq.num]] <- trimResults(res)
        }
        ##################
        #check whether all got candidates for each iteration
        ##################
        # number of clusters found in each of clusterings (candidates of length 2).
        # Note, number of clusters should be <= top.can (# top candidates to keep)
        nClusterPerK<-sapply(candidates,length) 
        if(any(nClusterPerK==0)){
            whyStop<-paste("Stopped in midst of searching for cluster",nfound+1," because no clusters meeting criteria found for iteration k=",k+1,"and previous clusters not similar enough.")
            if(verbose) cat(paste("Found ",paste(nClusterPerK,collapse=","),"clusters for k=",paste(k+seq_len(seq.num)-1,collapse=","),", respectively. Stopping iterating because no clusters meet trimming criteria found for this iteration and previous clusters not similar enough.\n"))
            break
        }
        if(isIntegerTop){ 
            whInvalid<-sweep(index.m,MARGIN=2,STATS=nClusterPerK, FUN=">")
            # vector (of length equal to number of combinations) saying which of the clusterings didn't have at least top.can clusters:
            whInvalid<-which(apply(whInvalid,1,any))
            if(length(whInvalid)==nrow(index.m)){
                #all invalid -- probably means that for some k there were no candidates found. So should stop.
                stop("Internal coding error -- shouldn't have reached this point in seqCluster")
            }
            if(length(whInvalid)>0){
                if(verbose) cat("Did not find", top.can,"clusters for all k: ")
                if(verbose) cat(sprintf("found %s clusters for k= %s, respectively\n",paste(nClusterPerK,collapse=","),paste(k+seq_len(seq.num)-1,collapse=",")))
                tempIndex<-index.m[-whInvalid,,drop=FALSE]
            }
            else tempIndex<-index.m
        }
        else{
            tempIndex<- as.matrix(expand.grid(lapply(nClusterPerK, 
                function(ii) seq_len(ii))))
        }

        ##################
        #Calculate the stability pairwise between all of cluster combinations
        ##################
        #function to calculate the stability across the given combination
        #	y is a combination (row of tempIndex) giving clusters to compare stability from k, k+1
        calc.beta <- function(y) {
            #written generally enough to deal with seq.num>2; could be a lot simpler with seq.num=2.
            temp <- lapply(seq_len(seq.num), function(z) candidates[[z]][[y[z]]]) #each
            i.temp <- temp[[1]]
            if(betaNum %in% c("all","first")) u.temp <- temp[[1]] ###eap: changed here
            if(betaNum == "last") u.temp<-temp[[seq.num]]
            for (j in 2:seq.num) {
                i.temp <- intersect(i.temp, temp[[j]])
                if(betaNum=="all") u.temp <- union(u.temp, temp[[j]])
            }
            out<-length(i.temp)/length(u.temp)
            return(out)
            
        }
        beta.temp <- apply(tempIndex, 1, calc.beta) 
        if (any(beta.temp >= beta)){
            found <- TRUE
            nfound <- nfound + 1
            if(verbose) cat(paste("Cluster",nfound,"found."), "")
            if (k.start > k.min) k.start <- k.start - 1 #decrease
            found.temp <- candidates[[whReturn]][[tempIndex[which.max(beta.temp)[1], whReturn]]]
            kend[[nfound]]<-k+seq.num-1 #just assuming returning last here!
            kstart[[nfound]]<-currentStart
            tclust[[nfound]]<-colnames(inputMatrix)[found.temp]
            tclust[[nfound]]<-as.numeric(tclust[[nfound]])
            if(input == "diss") inputMatrix<-inputMatrix[-found.temp, -found.temp, drop=FALSE]
            else inputMatrix<-inputMatrix[ , -found.temp, drop=FALSE]
            remain <- remain - length(tclust[[nfound]])
            if(verbose) cat(paste("Cluster size:", length(tclust[[nfound]]),
                                  "\tRemaining number of points:", remain, "\n"),
                            "")
        }
        else {
            found = FALSE
            k = k + 1
        }
    }
    ###############
    # Clean up and return results
    ###############
    if(is.null(whyStop)){
        if(remain< remain.n) whyStop<-"Ran out of samples"
        if(!found & k>k.max) whyStop<-paste("Went past k.max=",k.max,"in looking for cluster with similarity to previous.")
    }
    clusterVector<-.clusterListToVector(tclust,N)
    if(all(clusterVector==-1) & length(tclust)>0) stop("coding error")
    if(nfound>0){
        size <- sapply(tclust, length)
        sizeMat<-cbind(size=size,kStart=kstart,kEnd=kend,nIter=kend-kstart)
        res <- list(clustering = clusterVector, clusterInfo = sizeMat, whyStop=whyStop)
        if(verbose) cat(paste("Stopped because:", whyStop),"")
        return(res)
    }
    else{
        if(verbose) cat("No tight clusters could be found with given parameters")
        return(list(clustering = clusterVector, whyStop=whyStop))
    }
    
}
