#' @title Make hierarchy of set of clusters
#'
#' @aliases makeDendrogram,ClusterExperiment-method
#' @description Makes a dendrogram of a set of clusters based on hclust on the
#'   medoids of the cluster.
#' @param x data to define the medoids from. Matrix and
#'   \code{\link{ClusterExperiment}} supported.
#' @param cluster A numeric vector with cluster assignments. If x is a
#'   ClusterExperiment object, cluster is automatically the primaryCluster(x).
#'   ``-1'' indicates the sample was not assigned to a cluster.
#' @param unassignedSamples how to handle unassigned samples("-1") ; only
#'   relevant for sample clustering. See details.
#' @param reduceMethod character A character identifying what type of
#'   dimensionality reduction to perform before clustering. Can be either a
#'   value stored in either of reducedDims or filterStats slot or a built-in
#'   diminsionality reduction/filtering. The option "coCluster" will use the
#'   co-Clustering matrix stored in the 'coClustering' slot of the
#'   \code{ClusterExperiment} object
#' @param nDims The number of dimensions to keep from \code{reduceMethod}. If
#'   missing calls \code{\link{defaultNDims}}.
#' @param whichCluster an integer index or character string that identifies
#'   which cluster should be used to make the dendrogram. Default is
#'   primaryCluster.
#' @param ... for makeDendrogram, if signature \code{matrix}, arguments passed
#'   to hclust; if signature \code{ClusterExperiment} passed to the method for
#'   signature \code{matrix}. For plotDendrogram, passed to
#'   \code{\link{plot.dendrogram}}.
#' @inheritParams clusterSingle
#' @inheritParams reduceFunctions
#' @details The function returns two dendrograms (as a list if x is a matrix or
#'   in the appropriate slots if x is ClusterExperiment). The cluster dendrogram
#'   is created by applying \code{\link{hclust}} to the medoids of each cluster.
#'   In the sample dendrogram the clusters are again clustered, but now the
#'   samples are also part of the resulting dendrogram. This is done by giving
#'   each sample the value of the medoid of its cluster.
#' @details The argument \code{unassignedSamples} governs what is done with
#'   unassigned samples (defined by a -1 cluster value). If
#'   unassigned=="cluster", then the dendrogram is created by hclust of the
#'   expanded medoid data plus the original unclustered observations. If
#'   \code{unassignedSamples} is "outgroup", then all unassigned samples are put
#'   as an outgroup. If the \code{x} object is a matrix, then
#'   \code{unassignedSamples} can also be "remove", to indicate that samples
#'   with "-1" should be discarded. This is not a permitted option, however,
#'   when \code{x} is a \code{ClusterExperiment} object, because it would return
#'   a dendrogram with less samples than \code{NCOL(x)}, which is not permitted
#'   for the \code{@dendro_samples} slot.
#' @details If any merge information is stored in the input object, it will be
#'   erased by a call to makeDendrogram.
#' @return If x is a matrix, a list with two dendrograms, one in which the
#'   leaves are clusters and one in which the leaves are samples. If x is a
#'   ClusterExperiment object, the dendrograms are saved in the appropriate
#'   slots.
#'
#' @export
#' @seealso makeFilterStats, makeReducedDims
#' @examples
#' data(simData)
#'
#' #create a clustering, for 8 clusters (truth was 3)
#' cl <- clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#'
#' #create dendrogram of clusters:
#' hcl <- makeDendrogram(cl)
#' plotDendrogram(hcl)
#' plotDendrogram(hcl, leafType="samples",plotType="colorblock")
#'
#' @name makeDendrogram
#' @rdname makeDendrogram
#' @importFrom phylobase tipLabels
setMethod(
    f = "makeDendrogram",
    signature = "ClusterExperiment",
    definition = function(x, whichCluster="primaryCluster",reduceMethod="mad",
                          nDims=defaultNDims(x,reduceMethod),filterIgnoresUnassigned=TRUE,
                          unassignedSamples=c("outgroup", "cluster"),
                          whichAssay=1,...)
    {
        passedArgs<-list(...)
        
        checkIgnore<-.depricateArgument(passedArgs=passedArgs,"filterIgnoresUnassigned","ignoreUnassignedVar")
        if(!is.null(checkIgnore)){
            passedArgs<-checkIgnore$passedArgs
            filterIgnoresUnassigned<-checkIgnore$val
        }
        unassignedSamples<-match.arg(unassignedSamples)
        whCl<-.convertSingleWhichCluster(x,whichCluster,passedArgs)
        cl<-clusterMatrix(x)[,whCl]
        ##erase merge information
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
		#Don't really need this any more...
        x@dendro_outbranch<- .hasOutBranch(x)
        ch<-.checkDendrogram(x)
        if(!is.logical(ch)) stop(ch)
        return(x)
})



#' @rdname makeDendrogram
#' @export
setMethod(
    f = "makeDendrogram",
    signature = "dist",
    definition = function(x, cluster,
                          unassignedSamples=c("outgroup", "cluster", "remove"),
                          calculateSample=TRUE,...) {
        unassigned <- match.arg(unassignedSamples)
        if(is.null(attributes(x)$Labels)) {
            attributes(x)$Labels <- .makeSampleNames(seq_len(nSamples))
        }
		if(!all(is.na(suppressWarnings(as.numeric(attributes(x)$Labels ))))){
			warning("Cannot use the attributes(x)$Labels because they are numbers. Making sample names")
			sampleNames<-.makeSampleNames(seq_len(nSamples))
		}
		
        clusterD<-.makeClusterDendro(x,cluster,type="dist",...)  
        fullD<-.makeSampleDendro(x,dendro=clusterD, cl=.convertToNum(cluster), type=c("dist"), unassignedSamples=unassigned,sampleEdgeLength=0,  outbranchLength=1,calculateSample=calculateSample)
        return(list(samples=fullD,clusters=clusterD))
    })

#' @rdname makeDendrogram
#' @importFrom DelayedArray DelayedArray
#' @export
setMethod(
    f = "makeDendrogram",
    signature = "matrixOrHDF5",
    definition = function(x, cluster,
                          unassignedSamples=c("outgroup", "cluster", "remove"),
                          calculateSample=TRUE,...) {
		
        unassigned <- match.arg(unassignedSamples)
        if(is.null(colnames(x))) {
            colnames(x) <- .makeSampleNames(seq_len(ncol(x)))
        }
		if(!all(is.na(suppressWarnings(as.numeric(colnames(x) ))))){
			warning("Cannot use the colnames(x) because they are numbers. Making sample names")
			sampleNames<-.makeSampleNames(seq_len(ncol(x)))
		}
        clusterD<-.makeClusterDendro(x,cluster,type="mat",...)
        fullD<- .makeSampleDendro(x,clusterDendro=clusterD, cl=.convertToNum(cluster), type=c("mat"), unassignedSamples=unassigned, sampleEdgeLength=0,  outbranchLength=1,calculateSample=calculateSample)
        
        return(list(samples=fullD,clusters=clusterD))
    })

#' @importFrom phylobase nodeLabels<- nodeLabels getNode phylo4d
.makeClusterDendro<-function(x, cl,type=c("mat","dist"),...){
    type<-match.arg(type)
    if(type=="dist"){
        nSamples<-  attributes(x)$Size
        if(nSamples != length(cl)) {
            stop("cluster must be the same length as the number of samples")
        }
    }
    if(type=="mat"){
        if(ncol(x) != length(cl)) {
            stop("cluster must be the same length as the number of samples")
        }
     }
    
    clNum<-.convertToNum(cl)
    
    #############
    # Cluster dendrogram
    #############
    whKeep <- which(clNum >= 0) #remove -1, -2
    if(length(whKeep) == 0) stop("all samples have clusterIds<0")
    if(length(unique(cl[whKeep]))==1) stop("Only 1 cluster given. Can not make a dendrogram.")
    clFactor <- factor(cl[whKeep])
    
    #each pair of clusters, need to get median of the distance values
    #do a double by, just returning the values as a vector, and then take the median
    if(type=="dist"){
        medoids<-do.call("rbind", by(as.matrix(x)[whKeep,whKeep], clFactor, function(z){
            out<-as.vector(by(t(z),clFactor,function(y){median(as.vector(unlist(y)))}))
            names(out)<-levels(clFactor)
            return(out)
        }))
        diag(medoids)<-0 #incase the median of within is not zero...
        rownames(medoids) <- levels(clFactor)
        colnames(medoids) <- levels(clFactor)
    }
    if(type=="mat"){
        medoids <- do.call("rbind", by(t(x[,whKeep]), clFactor, function(z){apply(z, 2, median)}))
        rownames(medoids) <- levels(clFactor)
    }
    nPerCluster <- table(clFactor)
    clusterD<-if(type=="dist").convertToPhyClasses(stats::hclust(as.dist(medoids),members=nPerCluster,...),returnClass=c("phylo4")) else .convertToPhyClasses(stats::hclust(dist(medoids)^2,members=nPerCluster,...),returnClass=c("phylo4"))
       
   #create data for phylo4d object
   clusterNodes<-phylobase::getNode(clusterD,type="all") 
   clusterIdDendro<-paste("ClusterId",names(clusterNodes),sep="")
   clusterIdDendro[is.na(names(clusterNodes))]<-NA
   
   #make permanent node ids, starting with internal nodes.
   justNodes<-phylobase::getNode(clusterD,type="internal")
   justTips<-phylobase::getNode(clusterD,type="tip")
   nodeId<-rep(NA,length(clusterNodes))  #want tips last, and internal nodes first
   nodeId[match(justNodes,clusterNodes)]<-1:length(justNodes)
   nodeId[match(justTips,clusterNodes)]<-1:length(justTips)+length(justNodes)
   if(!all(sort(nodeId)==sort(1:length(nodeId)))) stop("coding error in giving node names")
   
   data.cl <- data.frame(NodeId = paste("NodeId",nodeId,sep=""), ClusterIdDendro = clusterIdDendro, ClusterIdMerge= rep(NA,length(clusterNodes)),stringsAsFactors=FALSE)
	data.cl$Position<-factor(rep("cluster hierarchy node",nrow(data.cl)), levels=.positionLevels)
	row.names(data.cl)<-as.character(unname(clusterNodes))
    data.cl$Position[phylobase::getNode(clusterD,type="tip")]<-"cluster hierarchy tip"
	  
	#give default node labels:
	phylobase::labels(clusterD)[justNodes]<-data.cl$NodeId[justNodes]
	clusterD<-phylobase::phylo4d(x=clusterD, all.data = data.cl)


	# -- Position: one of either "cluster hierarchy node","cluster hierarchy tip","tip hierarchy","assigned tip","outbranch hierarchy node","unassigned tip","outbranch root"). This column should be internal and validity check that numbers correspond to clustering from @dendro_index. Needs to be a check on "ClusterExperiment", not the "clusterPhylo4d"
    # -- ClusterIdDendro (only for cluster): cluster Id in @dendro_index that corresponds to node (NA otherwise)
    # -- ClusterIdMerge  (only for cluster): cluster Id in @merge_index that corresponds to node (NA otherwise)
    # -- MatchToClusterHier (only for sample tree): for those that are part of the cluster hierarchy, the (permanent) node id from cluster hierarchy -- ie not the name but the number so can always link them. Use this to create checks that the node names are the same, grab the default names, etc. 
    # -- SampleIndex (only for sample tree): index to the columns in the assay
    
    
    return(clusterD)
}



#' @importFrom phylobase rootNode descendants tdata
.makeSampleDendro<-function(x,clusterDendro, cl,type=c("mat","dist"), unassignedSamples=c("remove","outgroup","cluster"),sampleEdgeLength=0, outbranchLength=0,calculateSample=TRUE){
    unassignedSamples<-match.arg(unassignedSamples)
    type<-match.arg(type)
    sampleNames<-if(type=="mat") colnames(x) else attributes(x)$Labels
	
    if(calculateSample){
        whPos<-which(cl>0) #this is copy close to length of n
        if(!is.null(sampleNames) && length(sampleNames)!=length(cl)) stop("sampleNames must be same length as cluster vector")
					#converts internal node id and cluster id to the node, tip labels respectively.
        phyloObj <- .convertToPhyClasses(clusterDendro,"phylo",convertNode=TRUE,convertTip=TRUE)
        if(!is.ultrametric(phyloObj)) stop("coding error -- the cluster dendrogram is not ultrametric")
        nSamples<-switch(type,"mat"=ncol(x),"dist"=attributes(x)$Size)
        
        ###########################
        ## I. Make outlier tree if needed
        ###########################
        outbranchHclust<-NULL
        whNeg<-which(cl<0) #technically should be able to do with whPos, but a pain
        if(length(whNeg)>0){
	        outNames<- if(!is.null(sampleNames)) sampleNames[whNeg] else paste("OutSample",whNeg)
        	if(unassignedSamples=="outgroup" | length(whPos)==0){
                #--------
                # 1. Make outlier tree 
                #--------
                if(length(whNeg) > 5 | length(whPos)==0){ 
                    outlierDat <- if(type=="mat") x[,whNeg,drop=FALSE] else as.matrix(x)[whNeg,whNeg,drop=FALSE]
                    outbranchHclust <- if(type=="mat") stats::hclust(dist(t(outlierDat))^2) else  stats::hclust(as.dist(outlierDat))
                    outTree<- .convertToPhyClasses(x=outbranchHclust,"phylo") 
                    if(length(outTree$tip.label)!=length(whNeg)) stop("coding error - given hclust doesn't have correct number of tips.")
                }
                else{ #construct tree with just root and tips:
                    outTree<-.makeFakeBinary(tipNames=outNames,rootEdgeLength=1,edgeLength=.1)
                    #have to make it ultrametric -- since arbitrary doesn't matter.
                    if(length(outNames)>1) outTree<-.force.ultrametric(outTree)
                }
            }
            if(unassignedSamples=="cluster"){
                if(type=="dist")stop("cannot use unassigned='cluster' if input is a distance matrix")
                #code from assignUnassigned
                clFactor <- factor(cl[-whNeg])
                medoids <- do.call("rbind", by(t(x[,-whNeg]), clFactor, function(z){apply(z, 2, median)}))
                rownames(medoids) <- levels(clFactor)
                classif<-.genericClassify(x=x[,whNeg],centers=medoids) 
                cl[whNeg]<-classif
                whPos<-which(cl>0)
                whNeg<-which(cl<0)
                if(length(whPos)!=length(cl)) stop("coding error -- predicting unassigned missed some")
                unassignedSamples<-"remove"
            }			
        }      
        
        ###################
        ## II. Make tree with main samples:
        ###################
        if(length(whPos)>0){
            #---------
            #make list of phylo trees for each cluster:
            #---------
            fakePhyloList <- tapply(sampleNames[whPos], as.factor(paste("ClusterId",cl[whPos],sep="")), .makeFakeBinary, simplify=FALSE)
            
            #need to reorder so in order of the tips of phyloObj
            m<-match(phyloObj$tip.label,names(fakePhyloList))
            if(any(is.na(m))) stop("coding error -- names in dendrogram of clusters don't match cluster ids")
            fakePhyloList<-fakePhyloList[m]
            newPhylo<-.addTreesToTips(mainTree=phyloObj,tipTrees=fakePhyloList)
        }
        else newPhylo<-outTree 
        
        
        ###################
        ## III. Add outlier samples if unassignedSamples=="outgroup"
        ###################
        if(unassignedSamples == "outgroup" & length(whNeg)>0 & length(whPos)>0){
            if(length(whNeg)>1){
                newPhylo<-.mergePhylo(tree1=newPhylo,tree2=outTree,mergeEdgeLength=outbranchLength)
            }
            else{
                newPhylo<-.mergePhylo(tree1=newPhylo, tree2=outNames, mergeEdgeLength=outbranchLength)
            }
        }
        
        ###################
        ## IV. Convert to return format
        ###################
        newPhylo<-.convertToPhyClasses(newPhylo,"phylo4")
				
		#create data for phylo4d object
		nodeLabs<-phylobase::getNode(newPhylo,type="all") 
		mClusterHier<-rep(NA,length(nodeLabs))
		position<-rep(NA,length(nodeLabs))
		whInternal<-grep("NodeId",names(nodeLabs))
		mClusterHier[whInternal]<-as.character(.matchToDendroData(inputValue=names(nodeLabs)[whInternal], dendro=clusterDendro, columnValue="NodeId"))
		position[whInternal]<-as.character(.matchToDendroData(inputValue=names(nodeLabs)[whInternal], dendro=clusterDendro, columnValue="Position"))

		whCluster<-grep("ClusterId",names(nodeLabs))
		mClusterHier[whCluster]<-as.character(.matchToDendroData(inputValue=names(nodeLabs)[whCluster], dendro=clusterDendro, matchValue="ClusterIdDendro",columnValue="NodeId"))
		position[whCluster]<-as.character(.matchToDendroData(inputValue=names(nodeLabs)[whCluster], dendro=clusterDendro, matchValue="ClusterIdDendro",columnValue="Position"))

				#rather than keep track of position as make trees (which require making it a phylo4d class much more frequently), figure it out from here.
				
		rootNode<-phylobase::rootNode(newPhylo)
		rootChildren <- phylobase::descendants(newPhylo,node=rootNode,type="children")
			 
		 whNAChild<-which(is.na(names(rootChildren)))
		 if(length(whNAChild)==0){
			 #all are part of cluster hierarchy
			 #all NA nodes -> "tip hierarchy"
			 #all tips -> "assigned tip"
			 position[is.na(names(nodeLabs))]<-"tip hierarchy"
			 position[phylobase::getNode(newPhylo,type="tip")]<-"assigned tip"
			 #check if there is a singleton sample for outlier:
			 whSingletonOutlier<-which(is.na(mClusterHier[rootChildren]))
			 if(length(whSingletonOutlier)>0){
			 	position[rootChildren[whSingletonOutlier]]<-"unassigned tip"
			 }
	 
		 }
		 else if(length(whNAChild)==1){
			 #all nodes descending from child of NA -> "outbranch hierarchy node"
			 #all tips descending from child of NA -> "unassigned tip"
			 #all tips descending from non-NA -> "assigned tip"
			 #all NA nodes descending from non-NA -> "tip hierarchy"
			 #root -> "outbranch root"
			 position[rootNode]<-"outbranch root"
			 naNode<-rootChildren[which(is.na(names(rootChildren)))]
			 position[phylobase::descendants(newPhylo,node=naNode,type="ALL")] <- "outbranch hierarchy node" #includes tips, but will overwrite this in next line...
			 position[phylobase::descendants(newPhylo,node=naNode,type="tip")] <- "unassigned tip"
			 nonNANode<-rootChildren[which(!is.na(names(rootChildren)))]
			  whDescendantNonNA <- phylobase::descendants(newPhylo,node=nonNANode,type="all")
				whDescendantNonNA<-whDescendantNonNA[is.na(names(whDescendantNonNA))]
				position[whDescendantNonNA]<-"tip hierarchy"
			  whTipsDescendantNonNA <- phylobase::descendants(newPhylo,node=nonNANode,type="tip")
				position[whTipsDescendantNonNA]<-"assigned tip"
		
		 }
		else stop("coding error -- neither child of root is part of cluster hierarchy")
		if(any(is.na(position))) stop("coding error -- not all nodes assigned a position value")
		position<-factor(position,levels=.positionLevels)

		mSamples<-match(names(nodeLabs),sampleNames)
		#this will miss singleton clusters (which will have name ClusterIdX instead of sample name)
		whMissed<-intersect(which(is.na(mSamples)),grep("assigned tip",position))
		if(length(whMissed)>0){
			clusters<-.matchToDendroData(inputValue=mClusterHier[whMissed],clusterDendro,columnValue="ClusterIdDendro",matchValue="NodeId")
			if(any(is.na(clusters))) stop("coding error -- didn't find singleton cluster")
			whSamplesSingle<-match(gsub("ClusterId","",clusters),as.character(cl))
			mSamples[whMissed]<-whSamplesSingle
		}
		if(any(is.na(mSamples[position %in% c("unassigned tip","assigned tip")]))) stop("coding error -- finding sample index failed")
		data.cl <- data.frame(NodeId=mClusterHier, Position=position, SampleIndex=mSamples, stringsAsFactors=FALSE)				
		row.names(data.cl)<-as.character(nodeLabs)

		# #Need to fix back up the nodeLabels so match that of clusterDendro (in conversion made them the internal names)
		# whMatchCluster<-which(!is.na(data.cl$NodeId))
		# mToCluster<-.matchToDendroData(inputValue=data.cl$NodeId[whMatchCluster],clusterDendro,columnValue="matchIndex",matchValue="NodeId")
		# phylobase::labels(newPhylo,type="all")[whMatchCluster]<-phylobase::labels(clusterDendro,type="all")[mToCluster]
		
		phylobase::labels(newPhylo,type="all")<-NA
		return(phylobase::phylo4d(x=newPhylo, all.data = data.cl))

    }
    else return(NULL)
}


