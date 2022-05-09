#' @title Make hierarchy of set of clusters
#' @description Extends the makeDendrogram function from clusterExperiment package to the output of RSEC.
#' @name makeDendrogram
#' @rdname makeDendrogram
#' @aliases makeDendrogram,ClusterExperiment-method
#' @export
#' @seealso makeFilterStats, makeReducedDims
setMethod(
    f = "makeDendrogram",
    signature = "RSECClass",
    definition = function(x, whichCluster="primaryCluster",reduceMethod="mad",
                          nDims=defaultNDims(x,reduceMethod),filterIgnoresUnassigned=TRUE,
                          unassignedSamples=c("outgroup", "cluster"),
                          whichAssay=1,...)
    {
        passedArgs<-list(...)
        
        checkIgnore<-.depricateArgument(passedArgs=passedArgs,"filterIgnoresUnassigned","ignoreUnassignedVar") #06/2018 added in BioC 3.8
        if(!is.null(checkIgnore)){
            passedArgs<-checkIgnore$passedArgs
            filterIgnoresUnassigned<-checkIgnore$val
        }
        unassignedSamples<-match.arg(unassignedSamples)
        whCl<-getSingleClusterIndex(x,whichCluster,passedArgs)
        cl<-clusterMatrix(x)[,whCl]
        if(!is.na(mergeClusterIndex(x)) || !is.na(x@merge_dendrocluster_index)) x<-.eraseMerge(x)
        
        ########
        ##Transform the data
        ########
        if(length(reduceMethod)>1) stop('makeDendrogram only takes one choice of "reduceMethod" as argument')
        if(reduceMethod!="coCluster"){
            #need to change name of reduceMethod to make it match the
            #clustering information if that option chosen.
            datList<-getReducedData(object=x,whichCluster=whCl,reduceMethod=reduceMethod,
                                    nDims=nDims,filterIgnoresUnassigned=TRUE,  whichAssay=whichAssay,returnValue="list")
            x<-datList$objectUpdate
            dat<-datList$dat
            
            outlist <- do.call("makeDendrogram",c(list(
                x=dat, 
                cluster=cl,calculateSample=TRUE,
                unassignedSamples=unassignedSamples),
                passedArgs))
        }
        else{
            if(is.null(x@coClustering)) stop("Cannot choose 'coCluster' if 'coClustering' slot is empty. Run makeConsensus before running 'makeDendrogram' or choose another option for 'reduceMethod'")
            if(is.null(dimnames(x@coClustering))) stop("This ClusterExperiment object was made with an old version of clusterExperiment and did not give dimnames to the coClustering slot.")
            outlist<-do.call("makeDendrogram",c(list(
                x=as.dist(1-x@coClustering),
                cluster=cl,calculateSample=TRUE,
                unassignedSamples=unassignedSamples),
                passedArgs)) 
        }
        
        #Add clusterNames as Ids to cluster and sample dendrogram.
		x@dendro_clusters <- outlist$clusters
		phylobase::tipLabels(x@dendro_clusters)<-NA #erase any labels of the tips, internal nodes already have the defaults.
		x@dendro_samples <- outlist$samples #labels should have been erased already
        x@dendro_index<-whCl
        ch<-.checkDendrogram(x)
        if(!is.logical(ch)) stop(ch)
        return(x)
})


