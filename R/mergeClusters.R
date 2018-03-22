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
#'   the value used in the \code{mergeClusters} command is plotted at each node. 
#'   If 'none', then no proportions will be added to the dendrogram, though the 
#'   dendrogram will be drawn. 'plotInfo' can also be one of
#'   the valid input to \code{mergeMethod} (even if that method is not the method chosen in
#'   \code{mergeMethod} argument). \code{plotInfo} can also show the information
#'   corresponding  to "adjP" with a fold-change cutoff, by giving a value to this 
#'   argument in  the form of "adjP_2.0", for example.
#' @param isCount logical as to whether input data is a count matrix. See 
#'   details.
#' @param plot logical as to whether to plot the dendrogram with the merge
#'   results
#' @param nodePropTable Only for matrix version. Matrix of results from previous
#'   run of \code{mergeClusters} as returned by matrix version of
#'   \code{mergeClusters}. Useful if just want to change the cutoff. Not
#'   generally intended for user but used internally by package.
#' @param calculateAll logical. Whether to calculate the estimates for all 
#'   methods. This reduces computation costs for any future calls to 
#'   \code{mergeClusters} since the results can be passed to future calls of
#'   \code{mergeClusters} (and for \code{ClusterExperiment} objects this is done
#'   automatically).
#' @param showWarnings logical. Whether to show warnings given by the methods. 
#'   The 'locfdr' method in particular frequently spits out warnings (which may 
#'   indicate that its estimates are not reliable). Setting 
#'   \code{showWarnings=FALSE} will suppress all warnings from all methods (not 
#'   just "locfdr"). By default this is set to \code{showWarnings=FALSE} by 
#'   default to avoid large number of warnings being produced by "locfdr", but
#'   users may want to be more careful to check the warnings for themselves.
#' @param logFCcutoff Relevant only if the \code{mergeMethod} selected is 
#' "adjP", in which case the calculation 
#' of the proportion of individual tests significant will also require that 
#' the estimated log-fold change of the features to be at least this large in 
#' absolute value. Value will be rounded to nearest tenth of an integer via
#' \code{round(logFCcutoff,digits=1)}. For any other method, this parameter is ignored.
#' @param ... for signature \code{matrix}, arguments passed to the 
#'   \code{\link{plot.phylo}} function of \code{ape} that plots the dendrogram. 
#'   For signature \code{ClusterExperiment} arguments passed to the method for 
#'   signature \code{matrix} and then if do not match those arguments, will be 
#'   passed onto \code{\link{plot.phylo}}.
#' @inheritParams clusterMany,matrix-method
#'   
#' @details \strong{Estimation of Proportion non-null} "Storey" refers to the
#'   method of Storey (2002). "PC" refers to the method of Pounds and Cheng
#'   (2004). "JC" refers to the method of Ji and Cai (2007), and implementation 
#'   of "JC" method is copied from code available on Jiashin Ji's website, 
#'   December 16, 2015 
#'   (http://www.stat.cmu.edu/~jiashun/Research/software/NullandProp/). "locfdr"
#'   refers to the method of Efron (2004) and is implemented in the package 
#'   \code{\link{locfdr}}. "MB" refers to the method of Meinshausen and Buhlmann
#'   (2005) and is implemented in the package \code{\link{howmany}}. "adjP" 
#'   refers to the proportion of genes that are found significant based on a FDR
#'   adjusted p-values (method "BH") and a cutoff of 0.05.
#' @details  \strong{Count correction} If  \code{isCount=TRUE}, and the input is
#'   a matrix, \code{log2(count + 1)} will be used for
#'   \code{\link{makeDendrogram}} and the original data with voom correction
#'   will be used in \code{\link{getBestFeatures}}). If input is
#'   \code{\link{ClusterExperiment}}, then setting \code{isCount=TRUE} also
#'   means that the log2(1+count) will be used as the transformation, like for
#'   the matrix case as well as the voom calculation, and will NOT use the 
#'   transformation stored in the object. If FALSE, then transformData(x) will be 
#'   given to the input and will be used for both \code{makeDendrogram} and 
#'   \code{getBestFeatures}, with no voom correction.
#'   
#' @details \strong{Control of Plotting} If \code{mergeMethod} is not equal to 
#'   'none' then the plotting will indicate where the clusters will be merged by
#'   making dotted lines of edges that are merged together (assuming 
#'   \code{plotInfo} is not 'none'). \code{plotInfo} controls simultaneously 
#'   what information will be plotted on the nodes as well as whether the dotted
#'   lines will be shown for the merged cluster. Notice that the choice of 
#'   \code{plotInfo} (as long as it is not 'none') has no effect on how the 
#'   dotted edges are drawn -- they are always drawn based on the 
#'   \code{mergeMethod}. If you choose \code{plotInfo} to not be equal to the 
#'   \code{mergeMethod}, then you will have a confusing picture where the dotted
#'   edges will be based on the clustering created by \code{mergeMethod} while
#'   the information on the nodes is based on a different method. Note that you
#'   can override \code{plotInfo} by setting \code{show.node.label=FALSE}
#'   (passed to plot.phylo), so that no information is plotted on the nodes, but
#'   the dotted edges are still drawn. If you just want plot of the dendrogram,
#'   with no merging performed nor demonstrated on the plot, see
#'   \code{\link{plotDendrogram}}.
#' @details If the dendrogram was made with option
#'   \code{unassignedSamples="cluster"} (i.e. unassigned were clustered in with
#'   other samples), then you cannot choose the option
#'   \code{leafType='samples'}. This is because the current code cannot reliably
#'   link up the internal nodes of the sample dendrogram to the internal nodes
#'   of the cluster dendrogram when the unassigned samples are intermixed. 
#' @return If `x` is a matrix, it returns (invisibly) a list with elements 
#'   \itemize{ 
#'		\item{\code{clustering}}{ a vector of length equal to ncol(x) 
#'   giving the integer-valued cluster ids for each sample. "-1" indicates the 
#'   sample was not clustered.} 
#'		\item{\code{oldClToNew}}{ A table of the old 
#'   cluster labels to the new cluster labels.} 
#'		\item{\code{nodeProp}}{ A table of
#'   the proportions that are DE on each node.} 
#'		\item{\code{nodeMerge}}{ A table of indicating for each node whether merged or not and the cluster id in the new clustering that corresponds to the node}
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
#' leafType="clusters",plotType="name")
#'
#' #compare merged to original
#' tableClusters(merged,whichClusters=c("mergeClusters","clusterSingle"))
#'
#' @export
#' @importFrom howmany howmany lowerbound
#' @importFrom locfdr locfdr
#' @rdname mergeClusters
setMethod(f = "mergeClusters",
          signature = signature(x = "matrix"),
          definition = function(x, cl, dendro=NULL, 
			  mergeMethod=c("none", "Storey","PC","adjP", "locfdr", "MB", "JC"),
			  plotInfo="none", 
			  nodePropTable=NULL, calculateAll=TRUE, showWarnings=FALSE,
              cutoff=0.1, plot=TRUE,isCount=TRUE, logFCcutoff=0, ...){  
  dendroSamples<-NULL #currently option is not implemented for matrix version...
  if(!is.numeric(logFCcutoff) || logFCcutoff<0) stop("Invalid value for the parameter 'logFCcutoff'")
  logFCcutoff<-round(logFCcutoff,digits=1)
  if(!is.numeric(cutoff) || cutoff<0 || cutoff > 1) stop("Invalid value for the parameter 'cutoff'.")
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
  ############
  #determine what methods asked to be calculated
  ############
  rmFCString<-function(x){
	  wh<-grep("adjP",x)
	  if(length(wh)>0){
		  x[wh]<-"adjP"
	  }
	  return(x)
  }
  addFCString<-function(fc){
	  paste("adjP", format(fc,nsmall=1),sep="_")
  }
  adjPFCMethod<-if(logFCcutoff>0) addFCString(logFCcutoff) else NULL
  if(calculateAll) whMethodCalculate=.availMergeMethods
  else{
	  whMethodCalculate<-if(!mergeMethod=="none") mergeMethod else c()
	  if(plotInfo=="all") whMethodCalculate<-.availMergeMethods
	  if(plotInfo%in% .availMergeMethods) whMethodCalculate<-unique(c(whMethodCalculate,plotInfo))  	
  }
  if(logFCcutoff>0 && "adjP" %in% whMethodCalculate){
	  whMethodCalculate<-c(whMethodCalculate,adjPFCMethod)
	  if(mergeMethod=="adjP") mergeMethod<-addFCString(logFCcutoff)

  }

  otherVals<-colnames(nodePropTable)[!colnames(nodePropTable)%in%c("Node","Contrast")]
  otherVals<-c(otherVals,adjPFCMethod)
  possibleValues<-unique(c("none", "all","mergeMethod",.availMergeMethods,otherVals))
  plotInfo<-match.arg(plotInfo,possibleValues)

  if(mergeMethod=="none" & plotInfo=="none" & !calculateAll) stop("mergeMethod and plotInfo both equal 'none'; nothing to be done.")
  if(plotInfo=="mergeMethod" & mergeMethod=="none") {
    stop("can only plot 'mergeMethod' results if one method is selected")
  }

  ############
  #determine whether need to calculate, or if already in nodePropTable
  ############
  ### check valid existing table
  if(!is.null(nodePropTable)){
    if(!all(c("Node","Contrast") %in% colnames(nodePropTable))) stop("nodePropTable must have columns with names 'Node' and 'Contrast'")
    if(!all(.availMergeMethods %in% colnames(nodePropTable))) stop("All of the methods' names must be included in colnames of nodePropTable (with NA if not calculated):", paste(.availMergeMethods,collapse=",",sep=""))
  }
  needCalculate<-is.null(nodePropTable) || any(!whMethodCalculate %in% names(nodePropTable)) || any(is.na(nodePropTable[,whMethodCalculate]))

  ############
  ### calculate the estimated proportions
  ############  
  if(needCalculate){
	  #get per-gene test-statistics for the contrasts corresponding to each node (and return all)
	  sigTable <- getBestFeatures(x, cl, 
		  contrastType=c("Dendro"), dendro=dendro,
		  contrastAdj=c("All"),
		  number=nrow(x), p.value=1, isCount=isCount)
	  #divide table into each node and calculate proportion.
	  sigByNode <- by(sigTable, sigTable$ContrastName, function(x) {
	      storey<-if("Storey" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_Storey,showWarnings=showWarnings) else NA
	      pc <-if("PC" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_PC,showWarnings=showWarnings) else NA
	      mb <-if("MB" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_MB,showWarnings=showWarnings) else NA
	      locfdr <-if("locfdr" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_locfdr,showWarnings=showWarnings) else NA
	      
		  jc <-if("JC" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_JC,showWarnings=showWarnings) else NA
		  
	      adjP<-if("adjP" %in% whMethodCalculate)  .m1_adjP(adjPvalues=x$adj,logFC=x$logFC, logFCcutoff=0) else NA
	      out<-c("Storey"=storey,"PC"=pc,"adjP"=adjP, "locfdr"=locfdr, "MB"=mb,"JC"=jc)
	      if(!is.null(adjPFCMethod) && adjPFCMethod %in% whMethodCalculate){
	      	adjPFC<-.m1_adjP(adjPvalues=x$adj,logFC=x$logFC, logFCcutoff=logFCcutoff)
			out<-c(out,adjPFC)
			names(out)[length(names(out))]<-adjPFCMethod
			return(out)
	      }
		  else return(out)
	  })
  }
  else{
  	.mynote(paste("Using existing results of per-node significance -- no new tests on individual features will be made."))
	  
    sigByNode<-by(nodePropTable,nodePropTable$Node,function(x){x})
  }
  newcl <- cl
  ############
  ### #go up tree and merge clusters
  ############  
  phylo4Obj <- .makePhylobaseTree(dendro, "dendro")
  if(mergeMethod != "none"){
    valsPerNode <- sapply(sigByNode, function(x) {signif(x[[mergeMethod]], 2)})
    nodesBelowCutoff <- names(valsPerNode)[which(valsPerNode<cutoff)] #names of nodes below cutoff
    
    #find nodes where *all* descendants are below cutoff
    allTipNames <- phylobase::labels(phylo4Obj)[phylobase::getNode(phylo4Obj, type=c("tip"))]
    whToMerge <- sapply(nodesBelowCutoff,function(node){
      desc <- phylobase::descendants(phylo4Obj, node, type = c("all"))
      return(all(names(desc) %in% nodesBelowCutoff | names(desc) %in% allTipNames))
    })
    if(length(whToMerge)>0 && length(which(whToMerge)) > 0){
      nodesToMerge <- nodesBelowCutoff[whToMerge]
      
      #now find top ones
      whAnc <- sapply(nodesToMerge, function(node){
        anc <- phylobase::ancestors(phylo4Obj, node, type="all")
        return(!any(names(anc) %in% nodesToMerge))
      })
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
      uc<-sort(unique(newcl[newcl>0]))
      newcl <- as.numeric(factor(newcl, levels=uc,labels=as.character(1:length(uc))))
      #deal with -1/-2
      newcl[is.na(newcl)] <- cl[is.na(newcl)]
    }
  }
  ############
  ### save/update merge information
  ############  
  nodePropTableGiven<-nodePropTable
  whNotProp <- which(c("Node", "Contrast") %in% names(nodePropTableGiven)) #interger 0 if is NULL
  if(!needCalculate){
    #------------	
	## Just split the table apart 
    #------------	
    nodePropTable <- nodePropTableGiven[, -whNotProp]
    annotTable <- nodePropTableGiven[, whNotProp]
  }
  else{
	#------------	
	#create new table with calculated results
	#------------	
	nodePropTable <- do.call("rbind", sigByNode)
	annotTable <- data.frame("Node"=names(sigByNode),
	                       "Contrast" = as.character(sigTable$Contrast[ match(names(sigByNode), sigTable$ContrastName) ]), 
						   stringsAsFactors =FALSE)
    if (!is.null(nodePropTableGiven)) {
	  #------------	
      #add results from nodePropTableGiven (original table) to nodePropTableGiven
	  #------------	
      #check same nodes and contrasts with original table
      if (!identical(unname(sort(nodePropTableGiven$Node)),unname(sort(annotTable$Node))))
        stop("different nodes in given `nodePropTable` than those just calculated")
      if (!identical(unname(sort(nodePropTableGiven$Contrast)), sort(unname(annotTable$Contrast))))
        stop("different contrast values in given `nodePropTable` than those just calculated")
	  #----------
	  # Find methods in nodePropTableGiven not calculated this time
	  #----------
	  whTransfer<-colnames(nodePropTableGiven)[-whNotProp]
	  #1. do those that shared between two tables, but all NA in current, and were calculated previously
	  whTransferExist<-whTransfer[whTransfer %in% colnames(nodePropTable)]
	  whTransferExist<-whTransferExist[sapply(whTransferExist,function(x){ all(!is.na(nodePropTableGiven[,x])) })]
      whTransferExist<-whTransferExist[which(!whTransferExist %in% whMethodCalculate)]
      if(length(whTransferExist)>0){
        m<-match(annotTable$Node,nodePropTableGiven$Node)
		
        nodePropTable[,whTransferExist]<-data.matrix(nodePropTableGiven[m,whTransferExist,drop=FALSE])
		
      }
	  #2. Add those in original table (that were calculated), but don't yet exist in current nodePropTable
	  whTransfer<-whTransfer[!whTransfer %in% colnames(nodePropTable)]
	  if(length(whTransfer)>0){
		  nodePropTable<-cbind(nodePropTable,
				data.matrix(nodePropTableGiven[,whTransfer,drop=FALSE])	  	
				)
	  }
	  if(any(duplicated(colnames(nodePropTable)))) stop("Internal Coding error -- have multiple copies of a method values.")
	  #3. Go ahead and sort those that are the logFC so nice looking table (and make sure at the end of the table)
	  whSort<-which(!colnames(nodePropTable) %in% .availMergeMethods)
	  if(length(whSort)>0){
		  ord<-order(colnames(nodePropTable)[whSort])
		  nodePropTable<-cbind(nodePropTable[,.availMergeMethods], nodePropTable[,whSort[ord],drop=FALSE])	  
	  }
    }
    
  }						  
  ############
  #add merge information:
  #also determine whether node corresponds to a cluster in merge clusters
  ############
  if (mergeMethod != "none" &&
      length(whToMerge) > 0 && length(which(whToMerge)) > 0) {
	logicalMerge <- annotTable$Node %in% nodesToMerge
    #gives the names of original cluster ids
    corrspCluster <- sapply(annotTable$Node, function(node) {
      tips <- phylobase::descendants(phylo4Obj, node, type = c("tips")) #names of tips
      if (any(!names(tips) %in% as.character(cl))) {
        stop("coding error-- tips don't match values of cl")
      }
      m <- match(names(tips), cl) #gives index of first match of tips to cl
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
  nodeMergeTable<-data.frame(annotTable,
				"isMerged"=logicalMerge,
				"mergeClusterId"=corrspCluster,
				stringsAsFactors=FALSE)
  nodePropTable<-data.frame(annotTable,
	  			nodePropTable,
				stringsAsFactors=FALSE)
  
  if(mergeMethod=="none"){
    newcl<-NULL #was just the original and nothing changed, so don't return something that makes it look like theres a new clustering
    oldClToNew<-NULL
  }
  else{
    oldClToNew=table(Original=cl, New=newcl)
    #check node identification from above
    nmerge<-apply(oldClToNew,2,function(x){sum(x>0)})
    clustersThatMerge<-colnames(oldClToNew)[which(nmerge>1)]
    tableMergedClusters<- unname(sort(as.character(na.omit(nodeMergeTable$mergeClusterId))))
	actualMergedClusters<-sort(unname(clustersThatMerge))
	if(!identical(tableMergedClusters,actualMergedClusters)) stop("coding error -- wrong identification of merged clusters")
  }
  out<-list(clustering=newcl, oldClToNew=oldClToNew, cutoff=cutoff, logFCcutoff=logFCcutoff,
            nodeProp=nodePropTable, nodeMerge=nodeMergeTable,originalClusterDendro=dendro,mergeMethod=mergeMethod)
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
	combTable<-.nodeMergeInfo(nodeProp=out$nodeProp,nodeMerge=out$nodeMerge)
	
    if(!is.null(dendroSamples)) .plotDendro(dendroSamples,leafType="samples",mergeOutput=combTable,mergePlotType=plotInfo,mergeMethod=mergeMethod,cl=cl,plotType="name",outbranch=any(cl<0),...)
    else .plotDendro(dendro,leafType="clusters",mergeOutput=combTable,mergePlotType=plotInfo,mergeMethod=mergeMethod,cl=clMat,plotType="name",...)
    
  }
  invisible(out)
}
)

#' @rdname mergeClusters
#' @export
#' @param clusterLabel a string used to describe the type of clustering. By 
#'   default it is equal to "mergeClusters", to indicate that this clustering is
#'   the result of a call to mergeClusters (only if x is a ClusterExperiment object)
#' @param plotType if plotting, then whether leaves of dendrogram should be
#'   labeled by rectangular blocks of color ("colorblock")  or with the names of
#'   the leaves ("name") (only if x is a ClusterExperiment object). 
#' @param leafType if plotting, whether the leaves should be the clusters or the
#'   samples. Choosing 'samples' allows for visualization of how many samples 
#'   are in the merged clusters (only if x is a ClusterExperiment object), which
#'   is the main difference between choosing "clusters" and "samples",
#'   particularly if \code{plotType="colorblock"}
#' @details When the input is a \code{ClusterExperiment} object, the function
#'     attempts to update the merge information in that object. This is done by
#'     checking that the existing dendrogram stored in the object  (and run on
#'     the clustering stored in the slot \code{dendro_index}) is the same
#'     clustering that is stored in the slot \code{merge_dendrocluster_index}.
#'     For this reason, new calls to \code{\link{makeDendrogram}} will erase the merge
#'     information saved in the object.
#' @details If \code{mergeClusters} is run with \code{mergeMethod="none"}, the
#'   function may still calculate the proportions per node if \code{plotInfo} is
#'   not equal to "none" or \code{calculateAll=TRUE}. If the input object was a
#'   \code{ClusterExperiment} object, the resulting information will be still
#'   saved, though no new clustering was created; if there was not an existing
#'   merge method, the slot \code{merge_dendrocluster_index} will be updated.

setMethod(f = "mergeClusters",
          signature = signature(x = "ClusterExperiment"),
          definition = function(x, eraseOld=FALSE,isCount=FALSE, 
             mergeMethod="none",
			 plotInfo="all",
			 clusterLabel="mergeClusters",
             leafType=c("samples","clusters"),
			 plotType=c("colorblock","name","ids"),
			 plot=TRUE,
			 ...) {
	plotType<-match.arg(plotType)
	leafType<-match.arg(leafType)
	if(is.null(x@dendro_clusters)) {
	stop("`makeDendrogram` needs to be called before `mergeClusters`")
	}
	else{
	cl<-clusterMatrix(x)[,dendroClusterIndex(x)]
	.mynote(paste("Merging will be done on '",clusterLabels(x)[dendroClusterIndex(x)],"', with clustering index",dendroClusterIndex(x)))
	}
	if(isCount) .mynote("If `isCount=TRUE` the data will be transformed with voom() rather than with the transformation function in the slot `transformation`. This makes sense only for counts.")
	if(!x@dendro_outbranch){
		if(any(cl<0) & leafType=="samples"){
			warning("You cannot set 'leafType' to 'samples' in plotting mergeClusters unless the dendrogram was made with unassigned/missing (-1,-2) set to an outgroup (see makeDendrogram)")
			leafType<-"clusters"
		}
	}
	###Divide ... into mergeCluster arguments and plotting arguments
	###Remove mergeClusters arguments in ... so can pass arguments to  
	mergeArgs<-.methodFormals("mergeClusters","matrix") #list with names equal to formal arguments
	passedArgs<-list(...)
	if(any(!names(passedArgs) %in% names(mergeArgs))){
	   plotArgs<-passedArgs[which(!names(passedArgs) %in% names(mergeArgs))]
	}
	else plotArgs<-NULL

	  
	  
###Note, plot=FALSE, and then manually call .plotDendro afterwards to allow for passage of colors, etc.
	##########
	#check to make sure all match, same, etc. 
	##########
	if(!is.na(x@merge_index)){
		if(x@merge_dendrocluster_index==x@dendro_index) propTable<-x@merge_nodeProp
		else propTable<-NULL
	}
	else propTable<-NULL
  outlist <- mergeClusters(x=if(!isCount) transformData(x) else assay(x),
                           cl=cl, nodePropTable=propTable,
                           dendro=x@dendro_clusters, plotInfo=plotInfo,plot=FALSE,
                           isCount=isCount,mergeMethod=mergeMethod, ...)
  nodeMerge<-outlist$nodeMerge
  if(mergeMethod!="none"){#
	  #######################
	  ##add a new cluster, but only if there was a merging
	  #######################

    didMerge<-any(apply(outlist$oldClToNew,2,function(x){sum(x>0)>1}))
    if(!didMerge) .mynote("merging with these parameters did not result in any clusters being merged.")
    newObj <- ClusterExperiment(x, outlist$clustering,
                                transformation=transformation(x),
                                clusterTypes="mergeClusters", 
                                checkTransformAndAssay=FALSE)
	#-----							
	#add "m" to name of cluster
	#-----							
    newObj<-.addPrefixToClusterNames(newObj,prefix="m",whCluster=1)
    clusterLabels(newObj) <- clusterLabel
	#-----							
    ##Check if pipeline already ran previously and if so increase
	#-----							
    x<-.updateCurrentWorkflow(x,eraseOld,"mergeClusters")
	if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x)
    else retval<-.addBackSEInfo(newObj=newObj,oldObj=x)
	#-----							
    #Add merge info to merge slots
	#Do manually here, because need joint object to dendro_index stuff
	# (otherwise get validity errors)
	#Note will save merge info, even if no change from merging.
	#-----							
    retval@merge_nodeProp<-outlist$nodeProp
    retval@merge_index<-1
    retval@merge_method<-outlist$mergeMethod
    retval@merge_dendrocluster_index<-retval@dendro_index #update here because otherwise won't be right number.
    retval@merge_cutoff<-outlist$cutoff
	#----- 
	##The above can change the internal coding of the merge clusters (???Why???)
	##Need to update the nodeMerge to reflect this before save to object.
	##Do this by matching the outlist$clustering to the new clustering
	#----- 
	if(didMerge){ 
		
		if(all(is.na(nodeMerge$mergeClusterId))) stop("internal coding error -- merging done but no non-NA value in 'mergeClusterId' value") #just in case. Should be at least 1
		origOldToNew<-outlist$oldClToNew
		if(ncol(origOldToNew)>1){ #otherwise, only 1 cluster left, and will always have right number
			#
			currOldToNew<-tableClusters(retval,whichClusters=c(dendroClusterIndex(retval),mergeClusterIndex(retval)))
			#-----
			##Match merge Ids in nodeMerge to columns of origOldToNew
			#-----
			#get non NAs in nodeMerge
			whNotNAMerge<-which(!is.na(nodeMerge$mergeClusterId))
			idsInMergeTable<-nodeMerge$mergeClusterId[whNotNAMerge]
			#make origOldToNew only these ids in this order
			origOldToNew<-origOldToNew[,match(idsInMergeTable,colnames(origOldToNew)),drop=FALSE]
			#make new merge table in same order as these ids by match columns to each other
			mCols<-match(apply(origOldToNew,2,paste,collapse=","),apply(currOldToNew,2,paste,collapse=","))
			#currOldToNew<-currOldToNew[,mCols]
			nodeMerge$mergeClusterId[whNotNAMerge]<-as.numeric(colnames(currOldToNew)[mCols])
			
		}
		
	}
    retval@merge_nodeMerge<-nodeMerge	
	#------------
    ##Align the colors between mergeClusters and combineMany
	#------------
    tryAlign<-try(plotClusters(retval,resetColors = TRUE, whichClusters=c("mergeClusters","combineMany"),plot=FALSE),silent=TRUE)
	if(!inherits(tryAlign,"try-error")) 
		retval<-tryAlign
	else .mynote(paste("Unable to align mergeClusters with combineMany clusters because of the following error:",tryAlign))
    
  }
  else{ 
	##############
	#When no merging done, still need to save merge info so don't have to redo it.
	#But do not update the clustering, etc from above.
	##############
    retval<-x
    saveNodeTable<-is.na(x@merge_index)
    if(saveNodeTable && !is.na(x@merge_index) && x@merge_dendrocluster_index!=x@dendro_index) saveNodeTable<-FALSE
    if(saveNodeTable ){
		
      retval@merge_nodeProp <-outlist$nodeProp
    }
    if(is.na(retval@merge_index)) retval@merge_dendrocluster_index<-retval@dendro_index
    
  }
  ch<-.checkMerge(retval)
  if(!is.logical(ch)) stop(ch)  
  if(plot){
    dend<- switch(leafType, "samples"=retval@dendro_samples, "clusters"=retval@dendro_clusters)
    # leg<-clusterLegend(retval)[[retval@dendro_index]]
    #     cl<-switch(leafType,"samples"=clusterMatrix(retval)[,retval@dendro_index],"clusters"=NULL)
    if(leafType=="samples" & mergeMethod!="none" & plotType=="colorblock"){
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
    
    if(plotType=="id") leg[,"name"]<-leg[,"clusterIds"]
  	label<-switch(plotType,"name"="name","colorblock"="colorblock","ids"="name")
  	outbranch<-FALSE
  	if(leafType=="samples" & any(cl<0)) outbranch<-retval@dendro_outbranch
		#if(leafType=="samples" & any(cl<0)) outbranch<-TRUE

  	# outbranch<-any(clusterMatrix(retval)[,retval@dendro_index]<0)
  	# cl<-clusterMatrix(retval,whichCluster=retval@dendro_index)
  	# rownames(cl)<-colnames(retval)
  	# dend<-ifelse(leafType=="samples", retval@dendro_samples,retval@dendro_clusters)
  	if(!"legend" %in% names(plotArgs)) plotArgs$legend<-"none"
  	do.call(".plotDendro",c(list(dendro=dend,leafType=leafType,mergeOutput=.nodeMergeInfo(outlist$nodeProp,outlist$nodeMerge),mergePlotType=plotInfo,mergeMethod=mergeMethod,cl=cl,clusterLegendMat=leg,plotType=label,outbranch=outbranch,removeOutbranch=outbranch),plotArgs))
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
    if(!is.na(x@merge_dendrocluster_index)){
      if(!is.na(x@merge_index)){
        nodeMerge<-x@merge_nodeMerge
      }
      else{ #if run calculate all, can have the prop but no merge index
        nodeMerge<-cbind(
			x@merge_nodeProp[,c("Node","Contrast")],
			matrix(NA,nrow=nrow(x@merge_nodeProp),ncol=2)
			)
        colnames(nodeMerge)<-c("Node","Contrast","isMerged","mergeClusterId")
      }
	  return(.nodeMergeInfo(x@merge_nodeProp,nodeMerge))
    }

    else return(NULL)
  }
)

.nodeMergeInfo<-function(nodeProp,nodeMerge){
	if(!nrow(nodeProp)==nrow(nodeMerge)) stop("Coding error -- nodeProp and nodeMerge had different number of rows")
    out<-merge(x=nodeMerge,y=nodeProp,by=c("Node","Contrast"),all=TRUE)
	if(!nrow(out)==nrow(nodeProp)) stop("Coding error -- merging of nodeProp and nodeMerge indicates not same values of 'Node' and 'Contrast'")
    row.names(out)<-NULL
	return(out)
	
}
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



.myTryFunc<-function(FUN,showWarnings,...){
# If warn is zero (the default) warnings are stored until the top–level function returns. If 10 or fewer warnings were signalled they will be printed otherwise a message saying how many were signalled. An object called last.warning is created and can be printed through the function warnings. If warn is one, warnings are printed as they occur. If warn is two or larger all warnings are turned into errors.
	if(!showWarnings) suppressWarnings(x<-try(FUN(...),silent=TRUE))
	else x<-try(FUN(...),silent=TRUE)
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
.m1_adjP<-function(adjPvalues,logFC,logFCcutoff){
	if(length(adjPvalues)!=length(logFC)) stop("coding error -- adjPvalues and logFC must be of same length")
   sum(adjPvalues<=0.05 & abs(logFC)>=logFCcutoff)/length(adjPvalues)
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

#' @param by indicates whether output from \code{getMergeCorrespond} should be
#'   a vector/list with elements corresponding to merge cluster ids or elements
#'   corresponding to the original clustering ids. See return value for details.
#' @rdname mergeClusters
#' @return \code{getMergeCorrespond} returns the correspondence between the
#'   merged cluster and its originating cluster. If \code{by="original"} returns
#'   a named vector, where the names of the vector are the cluster ids of the
#'   originating cluster and the values of the vector are the cluster ids of the
#'   merged cluster. If \code{by="merge"} the results returned are organized by
#'   the merged clusters. This will generally be a list, with the names of the
#'   list equal to the clusterIds of the merge clusters and the entries the
#'   clusterIds of the originating clusters. However, if there was no merging
#'   done (so that the clusters are identical) the output will be a vector like
#'   with \code{by="original"}.
#' @aliases getMergeCorrespond
#' @export
setMethod(
  f = "getMergeCorrespond",
  signature = "ClusterExperiment",
  definition = function(x,by=c("merge","original")){
	by<-match.arg(by)
	if(is.na(x@merge_index)) stop("there is no merge clustering in this object")
    mCl<-clusterMatrix(x)[,x@merge_index]
    dCl<-clusterMatrix(x)[,x@merge_dendrocluster_index]
	mUnique<-as.character(unique(mCl[mCl>0]))
	dUnique<-as.character(unique(dCl[dCl>0]))
	nDCl<-length(dUnique)
	nMCl<-length(mUnique)
	notSubsetMessage<-"Current merge cluster is no longer a merge of the cluster in this -- indicated clustering that creates merge is not a partition of the merge cluster"
	if(nDCl<nMCl) stop(notSubsetMessage)
	if(nMCl==1){ #all d clusters are in the one cluster
		if(by=="merge"){
			out<-list(dUnique)
			names(out)<-mUnique
		}
		else{
			out<-lapply(dUnique,function(z){return(mUnique)})
		}
		return(out)

	}
    tab<-table(mCl,dCl) #match them up #might not work if only 1 cluster in one of these...
	if(is.null(dim(tab))) stop("coding error -- should get matrix of tabulations")
	margin<-switch(by,"merge"=1,"original"=2)
	dimvalue<-switch(by,"merge"=2,"original"=1)
	dimnms<-dimnames(tab)[[dimvalue]]
	nn<-apply(tab,margin,function(z){
        mm<-dimnms[z>0]
        if(by=="original" & length(mm)>1) stop(notSubsetMessage)
        else if (length(mm)==0) stop("coding error -- should have a corresponding merge cluster")
        else return(mm)
	})
	names(nn)<-dimnames(tab)[[margin]]
	return(nn)
})


