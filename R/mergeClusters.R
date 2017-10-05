.availMergeMethods<-c("adjP", "locfdr", "MB", "JC","PC","Storey")	
#' @title Merge clusters based on dendrogram
#'   
#' @description Takes an input of hierarchical clusterings of clusters and 
#'   returns estimates of number of proportion of non-null and merges those 
#'   below a certain cutoff.
#'   
#' @aliases mergeClusters
#'   
#' @param x data to perform the test on. It can be a matrix or a 
#'   \code{\link{ClusterExperiment}}.
#' @param cl A numeric vector with cluster assignments to compare to. ``-1'' 
#'   indicates the sample was not assigned to a cluster.
#' @param dendro dendrogram providing hierarchical clustering of clusters in cl.
#'   If x is a matrix, then the default is \code{dendro=NULL} and the function 
#'   will calculate the dendrogram with the given (x, cl) pair using 
#'   \code{\link{makeDendrogram}}. If x is a \code{\link{ClusterExperiment}} 
#'   object, the dendrogram in the slot \code{dendro_clusters} will be used. In 
#'   this case, this means that \code{\link{makeDendrogram}} needs to be called 
#'   before \code{mergeClusters}.
#' @param mergeMethod method for calculating proportion of non-null that will be
#'   used to merge clusters (if 'none', no merging will be done). See details 
#'   for description of methods.
#' @param cutoff minimimum value required for NOT merging a cluster, i.e. two 
#'   clusters with the proportion of DE below cutoff will be merged. Must be a 
#'   value between 0, 1, where lower values will make it harder to merge 
#'   clusters.
#' @param plotInfo what type of information about the merging will be shown on
#'   the dendrogram. If 'all', then all the estimates of proportion non-null
#'   will be plotted at each node of the dendrogram; if 'mergeMethod', then only
#'   the value used in the merging is plotted at each node. If 'none', then no
#'   proportions will be added to the dendrogram. 'plotInfo' can also be one of
#'   the mergeMethod choices (even if that method is not the method chosen in
#'   'mergeMethod' options).
#' @param isCount logical as to whether input data is a count matrix. See 
#'   details.
#' @param plot logical as to whether to plot the dendrogram with the merge
#'   results
#' @param nodePropTable Only for matrix version. Matrix of results from previous
#'   run of \code{mergeClusters} as returned by matrix version of
#'   \code{mergeClusters}. Useful if just want to change the cutoff. Not
#'   generally intended for user but used internally by package.
#' @param ... for signature \code{matrix}, arguments passed to the 
#'   \code{\link{plot.phylo}} function of \code{ape} that plots the dendrogram. 
#'   For signature \code{ClusterExperiment} arguments passed to the method for 
#'   signature \code{matrix} and then onto \code{\link{plot.phylo}}.
#' @inheritParams clusterMany,matrix-method
#'   
#' @details If  \code{isCount=TRUE}, and the input is a matrix, \code{log2(count
#'   + 1)} will be used for \code{\link{makeDendrogram}} and the original data
#'   with voom correction will be used in \code{\link{getBestFeatures}}). If
#'   input is \code{\link{ClusterExperiment}}, then setting \code{isCount=TRUE}
#'   also means that the log2(1+count) will be used as the transformation, like
#'   for the matrix case as well as the voom calculation, and will NOT use the 
#'   transformation stored in the object. If FALSE, then transform(x) will be 
#'   given to the input and will be used for both \code{makeDendrogram} and 
#'   \code{getBestFeatures}, with no voom correction.
#' @details "Storey" refers to the method of Storey (2002). "PC" refers to the 
#' method of Pounds and Cheng (2004). "JC" refers to the method of 
#' Ji and Cai (2007), and implementation 
#'   of "JC" method is copied from code available on Jiashin Ji's website, 
#'   December 16, 2015 
#'   (http://www.stat.cmu.edu/~jiashun/Research/software/NullandProp/). "locfdr"
#'   refers to the method of Efron (2004) and is implemented in the package 
#'   \code{\link{locfdr}}. "MB" refers to the method of Meinshausen and Buhlmann
#'   (2005) and is implemented in the package \code{\link{howmany}}. "adjP" 
#'   refers to the proportion of genes that are found significant based on a FDR
#'   adjusted p-values (method "BH") and a cutoff of 0.05.
#'   
#' @details If \code{mergeMethod} is not equal to 'none' then the plotting will 
#'   indicate where the clusters will be merged (assuming \code{plotInfo} is not
#'   'none'). Note setting both 'mergeMethod' and 'plotInfo' to 'none' will
#'   cause function to stop, because nothing is asked to be done. If you just
#'   want plot of the dendrogram, with no merging performed or demonstrated on
#'   the plot, see \code{\link{plotDendrogram}}.
#' @details If the dendrogram was made with option
#'   \code{unassignedSamples="cluster"} (i.e. unassigned were clustered in with
#'   other samples), then you cannot choose the option
#'   \code{leafType='samples'}. This is because the current code cannot reliably
#'   link up the internal nodes of the sample dendrogram to the internal nodes
#'   of the cluster dendrogram when the unassigned samples are intermixed. 
#' @return If `x` is a matrix, it returns (invisibly) a list with elements 
#'   \itemize{ \item{\code{clustering}}{ a vector of length equal to ncol(x) 
#'   giving the integer-valued cluster ids for each sample. "-1" indicates the 
#'   sample was not clustered.} \item{\code{oldClToNew}}{ A table of the old 
#'   cluster labels to the new cluster labels.} \item{\code{propDE}}{ A table of
#'   the proportions that are DE on each node.} 
#'   \item{\code{originalClusterDendro}}{ The dendrogram on which the merging 
#'   was based (based on the original clustering).} 
#'   \item{\code{cutoff}}{ The cutoff value for merging.} }
#' @return If `x` is a \code{\link{ClusterExperiment}}, it returns a new 
#'   \code{ClusterExperiment} object with an additional clustering based on the 
#'   merging. This becomes the new primary clustering.
#' @references Ji and Cai (2007), "Estimating the Null and the Proportion 
#' of Nonnull Effects in Large-Scale Multiple Comparisons", JASA 102: 495-906.
#' @references Efron (2004) "Large-scale simultaneous hypothesis testing: 
#' the choice of a null hypothesis," JASA, 99: 96-104.
#' @references Meinshausen and Buhlmann (2005) "Lower bounds for the 
#' number of false null hypotheses for multiple testing of associations", 
#' Biometrika 92(4): 893-907.
#' @references Storey (2002) "A direct approach to false discovery rates", J. R. Statist. Soc. B 64 (3)": 479-498.
#' @references Pounds and Cheng (2004). "Improving false discovery rate estimation." Bioinformatics 20(11): 1737-1745.

#' @seealso makeDendrogram, plotDendrogram, getBestFeatures
#' @examples
#' data(simData)
#' 
#' #create a clustering, for 8 clusters (truth was 3)
#' cl<-clusterSingle(simData, subsample=FALSE,
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#' 
#' #give more interesting names to clusters:
#' newNames<- paste("Cluster",clusterLegend(cl)[[1]][,"name"],sep="")
#' clusterLegend(cl)[[1]][,"name"]<-newNames
#' #make dendrogram
#' cl <- makeDendrogram(cl)
#'
#' #plot showing the before and after clustering
#' #(Note argument 'use.edge.length' can improve
#' #readability)
#' merged <- mergeClusters(cl, plotInfo="all",
#' mergeMethod="adjP", use.edge.length=FALSE)
#'
#' #Simpler plot with just dendrogram and single method
#' merged <- mergeClusters(cl, plotInfo="mergeMethod",
#' mergeMethod="adjP", use.edge.length=FALSE,
#' leafType="clusters",label="name")
#'
#' #compare merged to original
#' table(primaryCluster(cl), primaryCluster(merged))
#'
#' @export
#' @importFrom howmany howmany lowerbound
#' @importFrom locfdr locfdr
#' @rdname mergeClusters
setMethod(f = "mergeClusters",
          signature = signature(x = "matrix"),
          definition = function(x, cl, dendro=NULL,
                          mergeMethod=c("none", "Storey","PC","adjP", "locfdr", "MB", "JC"),
                          plotInfo=c("none", "all", "Storey","PC","adjP", "locfdr", "MB", "JC","mergeMethod"), nodePropTable=NULL,
                          cutoff=0.1, plot=TRUE,
                          isCount=TRUE,  ...) {  
  dendroSamples<-NULL #currently option is not implemented for matrix version...
  if(is.factor(cl)){
    warning("cl is a factor. Converting to numeric, which may not result in valid conversion")
    cl <- .convertToNum(cl)
  }
  if(!is.null(dendro)){
    #check valid
    ncluster <- length(table(cl[cl>0]))
    if(nobs(dendro) != ncluster) {
      stop("Not a valid input dendrogram (not equal to the number of non -1 clusters in cl).")
    }
  }
  mergeMethod <- match.arg(mergeMethod)
  plotInfo <- match.arg(plotInfo)
  if(mergeMethod=="none" & plotInfo=="none") stop("mergeMethod and plotInfo both equal 'none'; nothing to be done.")
  if(plotInfo=="mergeMethod" & mergeMethod=="none") {
    stop("can only plot 'mergeMethod' results if one method is selected")
  }
  #determine what methods asked to be calculated
  whMethodCalculate<-if(!mergeMethod=="none") mergeMethod else c()
  if(plotInfo=="all") whMethodCalculate<-.availMergeMethods
  if(plotInfo%in% .availMergeMethods) whMethodCalculate<-unique(c(whMethodCalculate,plotInfo))

	  #determine whether need to calculate, or if already in nodePropTable
  if(!is.null(nodePropTable)){
    if(!all(c("Node","Contrast") %in% colnames(nodePropTable))) stop("nodePropTable must have columns with names 'Node' and 'Contrast'")
    if(!all(.availMergeMethods %in% colnames(nodePropTable))) stop("All of the methods' names must be included in colnames of nodePropTable (with NA if not calculated):", paste(.availMergeMethods,collapse=",",sep=""))
  }
  needCalculate<-is.null(nodePropTable) || any(!whMethodCalculate %in% names(nodePropTable)) || any(is.na(nodePropTable[,whMethodCalculate]))
  if(needCalculate){### calculate the estimated proportions
	 
	  #get per-gene test-statistics for the contrasts corresponding to each node (and return all)
	  sigTable <- getBestFeatures(x, cl, contrastType=c("Dendro"), dendro=dendro,
	                               contrastAdj=c("All"),
	                              number=nrow(x), p.value=1, isCount=isCount)
	  #divide table into each node and calculate.
	  sigByNode <- by(sigTable, sigTable$ContrastName, function(x) {
	      storey<-if("Storey" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_Storey) else NA
	      pc <-if("PC" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_PC) else NA
	      mb <-if("MB" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_MB) else NA
	      locfdr <-if("locfdr" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_locfdr) else NA
	      jc <-if("JC" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_JC) else NA
	      adjP<-if("adjP" %in% whMethodCalculate)  .m1_adjP(x$adj) else NA
	      return(c("Storey"=storey,"PC"=pc,"adjP"=adjP, "locfdr"=locfdr, "MB"=mb,"JC"=jc))
	  })
  }
  else{
    sigByNode<-by(nodePropTable,nodePropTable$Node,function(x){x})
  }
  newcl <- cl
  phylo4Obj <- .makePhylobaseTree(dendro, "dendro")
  if(mergeMethod != "none"){
    #go up tree and merge clusters
    valsPerNode <- sapply(sigByNode, function(x) {signif(x[[mergeMethod]], 2)})
    nodesBelowCutoff <- names(valsPerNode)[which(valsPerNode<cutoff)] #names of nodes below cutoff

    #find nodes where *all* descendants are below cutoff
    allTipNames <- phylobase::labels(phylo4Obj)[phylobase::getNode(phylo4Obj, type=c("tip"))]
    whToMerge <- sapply(nodesBelowCutoff,function(node){
      desc <- phylobase::descendants(phylo4Obj, node, type = c("all"))
      return(all(names(desc) %in% nodesBelowCutoff | names(desc) %in% allTipNames))
    })
    #browser()
    if(length(whToMerge)>0 && length(which(whToMerge)) > 0){
      nodesToMerge <- nodesBelowCutoff[whToMerge]

      #now find top ones
      whAnc <- sapply(nodesToMerge, function(node){
        anc <- phylobase::ancestors(phylo4Obj, node, type="all")
        return(!any(names(anc) %in% nodesToMerge))
      })
     # browser()
      nodesAtTop <- nodesToMerge[whAnc]

      #make new clusters
      temp <- lapply(nodesAtTop, function(node){
        tips <- phylobase::descendants(phylo4Obj, node, type="tips")
        if(any(!names(tips) %in% as.character(cl))) {
          stop("coding error-- tips don't match values of cl")
        }
        newcl[cl%in% tips] <<- as.numeric(tips[[1]])
      })
      #make consecutive integers
      newcl <- as.numeric(factor(newcl, levels=unique(cl[cl>0])))
      #deal with -1/-2
      newcl[is.na(newcl)] <- cl[is.na(newcl)]
    }
  }
  nodePropTableGiven<-nodePropTable
	if(!needCalculate){
	  whProp <- which(names(nodePropTable) %in% .availMergeMethods)
	  nodePropTable <- nodePropTableGiven[, whProp]
	  annotTable <- nodePropTableGiven[, c("Node", "Contrast")]
	}
  else{
    if (!is.null(nodePropTable)) {
      #add results to given nodePropTable
      nodePropTable <- do.call("rbind", sigByNode)
      annotTable <- data.frame(
        "Node" = names(sigByNode),
        "Contrast" = as.character(sigTable$Contrast[match(names(sigByNode), sigTable$ContrastName)]),
        stringsAsFactors = FALSE
      )
      #check same nodes and contrasts
      if (!all(sort(nodePropTableGiven$Node) == sort(annotTable$Node)))
        stop("different nodes in nodePropTable than when calculated fresh")
      if (!all(sort(nodePropTableGiven$Contrast) == sort(annotTable$Contrast)))
        stop("different contrast values in nodePropTable than when calculated fresh")
      whInGiven<-.availMergeMethods[sapply(.availMergeMethods,function(x){ all(!is.na(nodePropTableGiven[,x])) })]
      whInGiven<-whInGiven[which(!whInGiven %in% whMethodCalculate)]
      if(length(whInGiven)>0){
        m<-match(annotTable$Node,nodePropTableGiven$Node)
        nodePropTable[,whInGiven]<-nodePropTableGiven[m,whInGiven]
      }
    }
    else{
      nodePropTable <- do.call("rbind", sigByNode)
      annotTable <- data.frame("Node"=names(sigByNode),
                               "Contrast"=as.character(sigTable$Contrast[match(names(sigByNode), sigTable$ContrastName)]),stringsAsFactors =FALSE)
    }
  }						  
  #add merge information:
  #also determine whether node corresponds to a cluster in merge clusters
  if (mergeMethod != "none" &&
      length(whToMerge) > 0 && length(which(whToMerge)) > 0) {
    logicalMerge <- annotTable$Node %in% nodesToMerge
    corrspCluster <- sapply(annotTable$Node, function(node) {
      tips <- phylobase::descendants(phylo4Obj, node, type = c("tips"))
      if (any(!names(tips) %in% as.character(cl))) {
        stop("coding error-- tips don't match values of cl")
      }
      m <- match(tips, cl)
      if (length(unique(newcl[m])) == 1) return(unique(newcl[m]))
      else NA
    })
    #Need to decide if any of these are nested inside each other!
    #Should be the largest is the parent that was merged
    if(length(na.omit(corrspCluster))!=length(unique(na.omit(corrspCluster)))){
      uniqueCorr<-unique(na.omit(corrspCluster))
      correctNode<-sapply(uniqueCorr,function(x){
        nodes<-annotTable$Node[which(corrspCluster==x)]
        if(length(nodes)>1){
          ntips<-sapply(nodes,function(node){length(phylobase::descendants(phylo4Obj, node, type = c("tips")))})
          maxnode<-nodes[which.max(ntips)]
          #check true assumption, no weird cases
          maxdesc<-phylobase::descendants(phylo4Obj, maxnode, type = c("all"))
          if(!all(nodes[-which.max(ntips)] %in% names(maxdesc))) stop("coding error -- largest samples wasn't parent node")
          return(maxnode)
        }
        else return(nodes)
      })
      corrspCluster<-rep(NA,length(corrspCluster))
      corrspCluster[match(correctNode,annotTable$Node)]<-uniqueCorr
    }
  } else{
    logicalMerge <- rep(FALSE, length = nrow(annotTable))
    corrspCluster <- rep(NA, length = nrow(annotTable))
  }
  nodePropTable<-data.frame(annotTable,"isMerged"=logicalMerge,"mergeClusterId"=corrspCluster,nodePropTable,stringsAsFactors=FALSE)
  
  if(mergeMethod=="none"){
    newcl<-NULL #was just the original and nothing changed, so don't return something that makes it look like theres a new clustering
    oldClToNew<-NULL
  }
  else{
    oldClToNew=table(Original=cl, New=newcl)
    #check node identification from above
    nmerge<-apply(oldClToNew,2,function(x){sum(x>0)})
    clustersThatMerge<-colnames(oldClToNew)[which(nmerge>1)]
    if(!all(sort(as.character(na.omit(nodePropTable$mergeClusterId)))==sort(clustersThatMerge))) stop("coding error -- wrong identification of merged clusters")
  }
  out<-list(clustering=newcl, oldClToNew=oldClToNew, cutoff=cutoff,
            propDE=nodePropTable, originalClusterDendro=dendro,mergeMethod=mergeMethod)
  if(plot){
    clMat<-cbind(Original=cl, mergeCluster=newcl)
    if(!is.null(dendroSamples)){
      if(is.null(names(cl))){
        warning("dendroSamples argument will be ignored because cl does not have names to allow for linkage to the dendroSamples values")
        dendroSamples<-NULL
      } 
      else{
        rownames(clMat)<-names(cl)	  	
      }
    }
    if(!is.null(dendroSamples)) .plotDendro(dendroSamples,leafType="samples",mergeOutput=out,mergePlotType=plotInfo,mergeMethod=mergeMethod,cl=cl,label="name",outbranch=any(cl<0),...)
    else .plotDendro(dendro,leafType="clusters",mergeOutput=out,mergePlotType=plotInfo,mergeMethod=mergeMethod,cl=clMat,label="name",...)
    
  }
  invisible(out)
}
)


#' @rdname mergeClusters
#' @export
#' @param clusterLabel a string used to describe the type of clustering. By 
#'   default it is equal to "mergeClusters", to indicate that this clustering is
#'   the result of a call to mergeClusters (only if x is a ClusterExperiment object)
#' @param labelType if plotting, then whether leaves of dendrogram should be
#'   labeled by rectangular blocks of color ("colorblock")  or with the names of
#'   the leaves ("name") (only if x is a ClusterExperiment object). 
#' @param leafType if plotting, whether the leaves should be the clusters or the
#'   samples. Choosing 'samples' allows for visualization of how many samples 
#'   are in the merged clusters (only if x is a ClusterExperiment object), which
#'   is the main difference between choosing "clusters" and "samples",
#'   particularly if \code{labelType="colorblock"}
setMethod(f = "mergeClusters",
          signature = signature(x = "ClusterExperiment"),
          definition = function(x, eraseOld=FALSE,isCount=FALSE, 
                                mergeMethod="none",plotInfo="all",clusterLabel="mergeClusters",
                                leafType=c("samples","clusters" ),labelType=c("colorblock","name","ids"),
                                plot=TRUE,...) {
  labelType<-match.arg(labelType)
  leafType<-match.arg(leafType)
  if(is.null(x@dendro_clusters)) {
    stop("`makeDendrogram` needs to be called before `mergeClusters`")
  }
  else{
    cl<-clusterMatrix(x)[,dendroClusterIndex(x)]
    note("Merging will be done on '",clusterLabels(x)[dendroClusterIndex(x)],"', with clustering index",dendroClusterIndex(x))
  }
  if(isCount) note("If `isCount=TRUE` the data will be transformed with voom() rather than
with the transformation function in the slot `transformation`.
This makes sense only for counts.")
	if(!x@dendro_outbranch){
		if(any(cl<0) & leafType=="samples"){
			warning("You cannot set 'leafType' to 'samples' in plotting mergeClusters unless the dendrogram was made with unassigned/missing (-1,-2) set to an outgroup (see makeDendrogram)")
			leafType<-"clusters"
		}
	}
  
###Note, plot=FALSE, and then manually call .plotDendro afterwards to allow for passage of colors, etc.
  if(!is.na(x@merge_index)){
    #check to make sure all match, same, etc. 
    if(x@merge_dendrocluster_index==x@dendro_index) propTable<-x@merge_nodeProp
    else propTable<-NULL
  }
  else propTable<-NULL
  outlist <- mergeClusters(x=if(!isCount) transform(x) else assay(x),
                           cl=cl, nodePropTable=propTable,
                           dendro=x@dendro_clusters, plotInfo=plotInfo,plot=FALSE,
                           isCount=isCount,mergeMethod=mergeMethod, ...)
  propTable<-outlist$propDE[,c("Node","Contrast",.availMergeMethods)]
  mergeTable<-outlist$propDE[,c("Node","Contrast","isMerged","mergeClusterId")]
  ##Did anything change??
  
  if(mergeMethod!="none"){#only add a new cluster if there was a mergeMethod. otherwise, mergeClusters just returns original cluster!
    didMerge<-any(apply(outlist$oldClToNew,2,function(x){sum(x>0)>1}))
    if(!didMerge) note("merging with these parameters did not result in any clusters being merged.")
    newObj <- clusterExperiment(x, outlist$clustering,
                                transformation=transformation(x),
                                clusterTypes="mergeClusters", 
								checkTransformAndAssay=FALSE)
    #add "m" to name of cluster
    newObj<-.addPrefixToClusterNames(newObj,prefix="m",whCluster=1)
    clusterLabels(newObj) <- clusterLabel
    ##Check if pipeline already ran previously and if so increase
    x<-.updateCurrentWorkflow(x,eraseOld,"mergeClusters")
    if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x)
    else retval<-.addBackSEInfo(newObj=newObj,oldObj=x)
    #add merge slots manually here, because need joint object to dendro_index stuff, and other wise get validity errors
    retval@merge_nodeProp<-propTable
    retval@merge_index<-1
    retval@merge_method<-mergeMethod
    retval@merge_nodeMerge<-mergeTable
    retval@merge_dendrocluster_index<-retval@dendro_index #update here because otherwise won't be right number.
    retval@merge_cutoff<-outlist$cutoff
        ch<-.checkMerge(retval)
    if(!is.logical(ch) || !ch) stop(ch)
  }
  else{ #still save merge info so don't have to redo it.
    retval<-x
    retval@merge_nodeProp=outlist$propDE[,c("Node","Contrast",.availMergeMethods)]
  }
  if(plot){
    dend<- switch(leafType,"samples"=retval@dendro_samples,"clusters"=retval@dendro_clusters)
  	# leg<-clusterLegend(retval)[[retval@dendro_index]]
  	#     cl<-switch(leafType,"samples"=clusterMatrix(retval)[,retval@dendro_index],"clusters"=NULL)
	if(leafType=="samples" & mergeMethod!="none" & labelType=="colorblock"){
		whClusters<-c(retval@dendro_index,primaryClusterIndex(retval))
	  	leg<-clusterLegend(retval)[whClusters]
	    cl<-clusterMatrix(retval,whichClusters=whClusters)
		rownames(cl)<-if(!is.null(colnames(retval))) colnames(retval) else as.character(1:ncol(retval))
		
	}
	else{
	  	leg<-clusterLegend(retval)[[retval@dendro_index]]
	  	    cl<-switch(leafType,"samples"=clusterMatrix(retval)[,retval@dendro_index],"clusters"=NULL)
		if(leafType=="samples"){
			names(cl)<-if(!is.null(colnames(retval))) colnames(retval) else as.character(1:ncol(retval))
		}
		
	}

    if(labelType=="id") leg[,"name"]<-leg[,"clusterIds"]
  	label<-switch(labelType,"name"="name","colorblock"="colorblock","ids"="name")
  	outbranch<-FALSE
  	if(leafType=="samples" & any(cl<0)) outbranch<-retval@dendro_outbranch
		#if(leafType=="samples" & any(cl<0)) outbranch<-TRUE

  # outbranch<-any(clusterMatrix(retval)[,retval@dendro_index]<0)
  # cl<-clusterMatrix(retval,whichCluster=retval@dendro_index)
  # rownames(cl)<-colnames(retval)
  # dend<-ifelse(leafType=="samples", retval@dendro_samples,retval@dendro_clusters)
     .plotDendro(dendro=dend,leafType=leafType,mergeOutput=outlist,mergePlotType=plotInfo,mergeMethod=mergeMethod,cl=cl,clusterLegendMat=leg,label=label,outbranch=outbranch,removeOutbranch=outbranch,legend="none")
  }
  
  invisible(retval)
}
)


#' @rdname mergeClusters
#' @return \code{nodeMergeInfo} returns information collected about the nodes
#'   during merging as a data.frame with the following entries:
#' \itemize{ \item{\code{Node}}{ Name of the node} 
#' \item{\code{Contrast}}{The
#' contrast compared at each node, in terms of the cluster ids} 
#' \item{\code{isMerged}}{ Logical as to whether samples from that node which were
#' merged into one cluster during merging} 
#' \item{\code{mergeClusterId}}{ If a
#' node corresponds to a new, merged cluster, gives the cluster id it
#' corresponds to. Otherwise NA} 
#' \item{\code{...}}{The remaining columns give
#' the estimated proportion of genes differentially expressed for each method. A
#' column of NAs means that the method in question hasn't been calculated yet.}
#' }
#' @aliases nodeMergeInfo
#' @export
setMethod(
  f = "nodeMergeInfo",
  signature = "ClusterExperiment",
  definition = function(x) {
    if(!is.na(x@merge_index)){
      nodeProp<-x@merge_nodeProp
      nodeMerge<-x@merge_nodeMerge
      m<-match(nodeMerge$Node,nodeProp$Node)
      out<-(cbind(nodeProp[m,c("Node","Contrast")],nodeMerge[,c("isMerged","mergeClusterId")],nodeProp[m,.availMergeMethods]))
      row.names(out)<-NULL
      return(out)
      
    }
    else return(NULL)
  }
)

#' @rdname mergeClusters
#' @return \code{mergeCutoff} returns the cutoff used for the current merging.
#' @aliases mergeCutoff
#' @export
setMethod(
  f = "mergeCutoff",
  signature = "ClusterExperiment",
  definition = function(x) {
    x@merge_cutoff
  }
)
#' @rdname mergeClusters
#' @return \code{mergeMethod} returns the method used for the current merge.
#' @aliases mergeMethod
#' @export
setMethod(
  f = "mergeMethod",
  signature = "ClusterExperiment",
  definition = function(x) {
    x@merge_method
  }
)
#' @rdname mergeClusters
#' @return \code{mergeClusterIndex} returns the index of the clustering used for the current merge.
#' @aliases mergeClusterIndex
#' @export
setMethod(
  f = "mergeClusterIndex",
  signature = "ClusterExperiment",
  definition = function(x) {
    x@merge_index
  }
)



.myTryFunc<-function(FUN,...){
  x<-try(FUN(...),silent=TRUE)
  if(!inherits(x, "try-error")) return(x)
  else return(NA)
}

#functions for estimating m1/m, the proportion of **non-null**
.m1_Storey<-function(pvalues,lambda=0.5){
	m<-length(pvalues)
	num<-length(which(pvalues>lambda))
	pi0<-num/(1-lambda)/m
	return(1-pi0)

}
.m1_PC<-function(pvalues){
	pi0<-2*mean(pvalues)
	return(1-pi0)

}

.m1_MB<-function(pvalues){
  nCorrect<-max(howmany::lowerbound(howmany::howmany(pvalues))) #the most you can call correctly
  return(nCorrect/length(pvalues))
}
.m1_adjP<-function(adjP){
  sum(adjP<0.05)/length(adjP)
}
.m1_locfdr<-function(tstats){
  locfdrResults<-locfdr::locfdr(tstats,plot=0)#ignore issue of df of t-statistic -- topTable doesn't give it out, and with large samples won't matter.
  p0<-locfdrResults$fp0["mlest","p0"] #estimate proportion null; ignore estimate of variability for now
  return(1-p0)
}
.m1_JC<-function(tstats){
  #copied code from Jianshin's website
  musigma<-try(.EstNull.func(tstats),silent=TRUE)
  .epsest.func(tstats,musigma$mu,musigma$s) #gives proportion of non-null
}


.makeMergeDendrogram<-function(object){
	if(is.na(object@dendro_index)) stop("no dendrogram for this clusterExperiment Object")
  #should this instead just silently return the existing?
	if(is.na(object@merge_index)) stop("no merging was done for this clusterExperiment Object")
  if(object@merge_dendrocluster_index != object@dendro_index) stop("dendrogram of this object was made from different cluster than that of merge cluster.")
  whClusterNode<-which(!is.na(object@merge_nodeMerge[,"mergeClusterId"]))
  clusterNode<-object@merge_nodeMerge[whClusterNode,"Node"]
  clusterId<-object@merge_nodeMerge[whClusterNode,"mergeClusterId"]
  phylo4Obj <- .makePhylobaseTree(object@dendro_clusters, "dendro")
  newPhylo4<-phylo4Obj
  for(node in clusterNode){
    #first remove tips of children nodes so all children of node are tips
    desc<-phylobase::descendants(newPhylo4, node, type = c("all"))
    whDescNodes<-which(names(desc) %in% phylobase::nodeLabels(newPhylo4))
    if(length(whDescNodes)>0){
      tipNodeDesc<-unlist(phylobase::descendants(newPhylo4, desc[whDescNodes], type = c("tips")))
      newPhylo4<-phylobase::subset(newPhylo4,tips.exclude=tipNodeDesc,trim.internal =FALSE)

    }
    #redo to check fixed problem
    desc<-phylobase::descendants(newPhylo4, node, type = c("all"))
    whDescNodes<-which(names(desc) %in% phylobase::nodeLabels(newPhylo4))
    if(length(whDescNodes)>0) stop("coding error -- didn't get rid of children nodes...")
    #should only have tips now
    tipsRemove<-phylobase::descendants(newPhylo4, node, type = c("tips"))
    newPhylo4<-phylobase::subset(newPhylo4,tips.exclude=tipsRemove,trim.internal =FALSE)
  }
  #return(newPhylo4)
  #Now need to change tip name to be that of the merge cluster
  mCl<-clusterMatrix(object)[,object@merge_index]
  dCl<-clusterMatrix(object)[,object@dendro_index]
  tab<-table(mCl,dCl) #match them up #might not work if only 1 cluster in one of these...

  newTips<-currTips<-phylobase::tipLabels(newPhylo4)
  whOldCl<-which(currTips %in% colnames(tab))
  if(length(whOldCl)>0){
    subtab<-tab[,currTips[whOldCl],drop=FALSE]
    nn<-apply(subtab,2,function(x){
      mm<-rownames(subtab)[x>0];
      if(length(mm)>1){
        stop("coding error -- shouldn't have multiple matches")
      }
      else if (length(mm)==0) stop("coding error -- shouldn't have no matches")
      else return(mm)
      })
    newTips[whOldCl]<-nn
  }
  if(!all(clusterNode %in% currTips)) stop("coding error -- some cluster nodes didn't wind up as tips of new tree")
  mClusterNode<-match(clusterNode, currTips)
  newTips[mClusterNode]<-as.character(clusterId)
  names(newTips)<-names(currTips) #this is the internal numbering of the nodes
  phylobase::tipLabels(newPhylo4)<-newTips
  #just some checks didn't screw up
  whMergeNode<-which(currTips %in% clusterNode)
  if(length(intersect(whMergeNode,whOldCl))) stop("coding error -- should be no overlap bwettween merged node in tree tips and old clusters")
  if(length(union(whMergeNode,whOldCl))!= length(currTips)) stop("coding error -- all tips should be either old clusters of merged nodes")

  if(length(currTips)!= length(unique(mCl[mCl>0]))) stop("coding error -- number of tips of new tree not equal to the number of clusters in merged cluster")
return(newPhylo4)
  #convert back to dendrogram class and return
}
 
  
