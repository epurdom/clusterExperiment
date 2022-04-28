############
## Checks that make up validity of ClusterExperiment:
############
.checkAssays<-function(object){
    # ############
    # check assays
    # ############
    if(length(assays(object)) < 1) {
        return("There must be at least one assay slot.")
    }
    if(!is.numeric(assay(object))) {
        return("The data must be numeric.")
    }
    if(anyNA(assay(object))) {
        return("NA values are not allowed.")
    }
    return(TRUE)
}
.checkTransform<-function(object){
    #############
    #check transform
    #############
    tX <- try(transformData(object),silent=TRUE)
    if(inherits(tX, "try-error")){
        stop(paste("User-supplied `transformation` produces error on the input data
                 matrix:\n",tX))
    }
    if(anyNA(tX)) {
        return("NA values after transforming data matrix are not allowed.")
    }
    return(TRUE)
}

.checkClusterMatrix<-function(object){
    ############
    #Check clusterMatrix
    ############
    zeroRow<-NROW(object@clusterMatrix)==0 #could happen in subsetting
    if(!all(is.na((object@clusterMatrix))) &
       !(NROW(object@clusterMatrix) == NCOL(object))) {
        return("If present, `clusterMatrix` must have as many row as cells.")
    }
    if(!zeroRow & !is.numeric(object@clusterMatrix)) {
        return("`clusterMatrix` must be a numeric matrix.")
    }
    
    if(NCOL(object@clusterMatrix)!= length(object@clusterTypes)) {
        return("length of clusterTypes must be same as NCOL of the clusterMatrix")
    }
    
    if(NCOL(object@clusterMatrix)!= length(object@clusterInfo)) {
        return("length of clusterInfo must be same as NCOL of the clusterMatrix")
    }
    #check internally stored as integers
    testConsecIntFun<-function(x){
        whCl<-which(!x %in% c(-1,-2))
        uniqVals<-unique(x[whCl])
        return(all(sort(unname(uniqVals))==seq_along(uniqVals)))
    }
    if(!zeroRow){
        testConsecIntegers<-apply(object@clusterMatrix,2,testConsecIntFun)
        if(!all(testConsecIntegers)) return("the cluster ids in clusterMatrix must be stored internally as consecutive integer values")		
    }
    
    return(TRUE)
}
.checkClusterLabels<-function(object){
    if(!all(is.na(object@clusterMatrix))){
        if(is.null(colnames(object@clusterMatrix))) return("must have clusterLabels by assignment of column names for clusterMatrix")
        if(any(duplicated(colnames(object@clusterMatrix)))) return("cannot have duplicated clusterLabels (i.e. clusterMatrix must have unique column names)")
        
    }
    return(TRUE)
}

#these functions are checks where don't need the corresponding object information
.checkDendroClusterFormat<-function(dendro,checkLabels=TRUE){
    data.cl<-phylobase::tdata(dendro)
    if(!all(.clusterDendroColumns %in% names(data.cl) )){
        return("dendro_clusters must have data with column names:",paste(.clusterDendroColumns,sep=","))
    }
    tp<-phylobase::tipLabels(dendro)
    if(checkLabels && any(sort(as.numeric(gsub("T","",tp)))!=sort(seq_along(tp))))
        return("dendro_clusters cannot have labels for the tips; user-defined labels for the tips (i.e. clusters) should be stored in the clusterLegend")
    if(any(is.na(data.cl$Position))) 
        return("dendro_clusters cannot have NA values in Position variable")
    if( any(is.na(data.cl$NodeId))) 
        return("dendro_clusters cannot have NA values in Node Id variable")
    return(TRUE)
}
.checkDendroSamplesFormat<-function(dendro,checkLabels=TRUE){
    data.cl<-phylobase::tdata(dendro,type="all")
    all(names(data.cl)%in% .clusterSampleColumns)
    if(!all(.clusterSampleColumns %in% names(data.cl) )){
        return("dendro_samples must have data with column names:",paste(.clusterDendroColumns,sep=","))
    }
    if(any(is.na(data.cl$Position))) 
        return("dendro_samples cannot have NA values in Position variable")
    if(checkLabels && any(!is.na(phylobase::nodeLabels(dendro)))) 
        return("dendro_samples cannot have labels for the nodes; the labels for the nodes should be in the cluster dendrogram")
    if(checkLabels){
        tp<-phylobase::tipLabels(dendro)
        if(any(sort(as.numeric(gsub("T","",tp)))!=sort(seq_along(tp)))) return("dendro_samples cannot have labels for the tips; the labels for the tips should be the colnames of the object")
    } 
    data.cl<-phylobase::tdata(dendro,type="tip")
    if(any(is.na(data.cl$SampleIndex))) 
        return("dendro_samples cannot have NA values for tips in SampleIndex variable")
    return(TRUE)
    
}
.checkDendrogram<-function(object){
    ############
    ##Check dendrogram
    ############
    if(!is.null(object@dendro_clusters)){
        if(is.na(dendroClusterIndex(object))) return("if dendrogram slots are filled, must have corresponding dendro_index defined.")
        dcluster<-clusterMatrix(object)[,dendroClusterIndex(object)]
        if(phylobase::nTips(object@dendro_clusters) != max(dcluster)) {
            return("dendro_clusters must have the same number of leaves as the number of (non-negative) clusters")
        }
        ch<-.checkDendroClusterFormat(object@dendro_clusters)
        if(!is.logical(ch)) return(ch)
        
        #further checks that require full CE object:
        data.cl<-phylobase::tdata(object@dendro_clusters,type="all")
        if(any(!gsub("ClusterId","",na.omit(data.cl$ClusterIdDendro)) %in% as.character(object@clusterMatrix[,object@dendro_index]))) return("ClusterIdDendro information in dendrogram slot must match the corresponding cluster ids in clustering defined by dendro_index slot")
    }
    else{
        if(!is.null(object@dendro_samples)) return("dendro_clusters should not be null if dendro_samples is non-null") #if comment out now optional to have samples
    }
    if(!is.null(object@dendro_samples)){
        if(phylobase::nTips(object@dendro_samples) != NCOL(object)) {
            return("dendro_samples must have the same number of leaves as the number of samples")
        }
        ch<-.checkDendroSamplesFormat(object@dendro_samples)
        if(!is.logical(ch)) return(ch)
        
    }
    else{
        if(!is.null(object@dendro_clusters)) return("dendro_samples should not be null if dendro_clusters is non-null") #if commented out, makes it optional to have samples
    }
    
    return(TRUE)
}


.checkPrimaryIndex<-function(object){
    ############
    ## Check primary index
    ############
    if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
        #check primary index
        if(length(object@primaryIndex) != 1) {
            if(length(object@primaryIndex) == 0) return("If more than one set of clusterings, a primary cluster must be specified.")
            if(length(object@primaryIndex) > 0) return("Only a single primary index may be specified")
        }
        if(object@primaryIndex > NCOL(object@clusterMatrix) | object@primaryIndex < 1) {
            return("`primaryIndex` out of bounds.")
        }
    }
    return(TRUE)
}

.checkClusterTypes<-function(object){
    ############
    ## Check clusterTypes
    ############
    if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
        if(NCOL(object@clusterMatrix) != length(object@clusterTypes)) {
            return("`clusterTypes` must be the same length as NCOL of `clusterMatrix`.")
        }
        if(!is.null(names(object@clusterTypes))) return("clusterTypes should not have names")
    }
    return(TRUE)
}

#this check is just to check object, not whether matches clusters
.checkClusterLegendList<-function(clusterLegend,allowNames=TRUE,reqNames=c("clusterIds","color","name")){
    if(!is.list(clusterLegend)) return("clusterLegend must be a list")
    if(!allowNames & !is.null(names(clusterLegend))) return("clusterLegend should not have names")
    testIsMatrix <- sapply(clusterLegend, function(x) {!is.null(dim(x))})
    if(!all(testIsMatrix)) {
        return("Each element of `clusterLegend` list must be a matrix")
    }
    testNames<-sapply(clusterLegend,function(x){
        if(is.null(colnames(x))) return(FALSE)
        else{
            if(!all(reqNames %in% colnames(x))) return(FALSE)
            else return(TRUE)
        }
    })
    if(!all(testNames)) {
        return(paste("each element of `clusterLegend` must be matrix with column names defined, with at a minimum the names names:", paste(reqNames,collapse=",")))
    }
    testColorCols1 <- sapply(clusterLegend, function(x){is.character(x)})
    if(!all(testColorCols1)) {
        return("each element of `clusterLegend` must be matrix of character values")
    }
    return(TRUE)
}
#this check checks both object (calls .checkClusterLegend) and whether matches clusters
#make so can call on arbitrary clusterLegend...not need to be CE object
.checkClustersWithClusterLegend<-function(clusters,clusterLegend){
    #check structure of clusterLegend list -- for CE object, can't have names.
    legendCheck<-.checkClusterLegendList(clusterLegend,allowNames=FALSE,reqNames=c("clusterIds","color","name"))
    if(!is.logical(legendCheck)) return(legendCheck)
    
    #check matches clusters
    if(length(clusterLegend) != NCOL(clusters)) {
        return("`clusterLegend` must be list of same length as NCOL of
               `clusterMatrix`")
    }
    testColorRows <- sapply(clusterLegend, function(x){nrow(x)})
    testClusterMat <- apply(clusters, 2, function(x) {length(unique(x))})
    if(!all(testColorRows == testClusterMat)) {
        return("each element of `clusterLegend` must be matrix with number of rows equal to the number of clusters (including -1 or -2 values) in `clusterMatrix`")
    }
    testColorCols1 <- sapply(seq_along(clusterLegend), function(ii){
        col<-clusterLegend[[ii]]
        x<-clusters[,ii]
        y<-col[,"clusterIds"]
        if(is.numeric(x)){
            y<-as.numeric(y)
        }
        return(all(y %in% x) & all(x %in% y))
    })
    if( !all(testColorCols1)) {
        return("each element of `clusterLegend` must be matrix with column
             `clusterIds` matching the corresponding
             clusterMatrix values")
    }
    
    return(TRUE)
}
.checkClusterLegend<-function(object){
    if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
        return(.checkClustersWithClusterLegend(clusters=object@clusterMatrix,clusterLegend=object@clusterLegend))
    }
    return(TRUE)
}

.checkOrderSamples<-function(object){
    ####
    #test orderSamples
    ####
    if(length(object@orderSamples)!=NCOL(assay(object))) {
        return("`orderSamples` must be of same length as number of samples
	       (NCOL(assay(object)))")
    }
    if(any(!object@orderSamples %in% seq_len(NCOL(assay(object))))) {
        return("`orderSamples` must be values between 1 and the number of samples.")
    }
    return(TRUE)
}


