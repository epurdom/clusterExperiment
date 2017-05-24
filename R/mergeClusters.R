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
#' @param plotType what type of plotting of dendrogram. If 'all', then all the 
#'   estimates of proportion non-null will be plotted at each node of the
#'   dendrogram; if 'mergeMethod', then only the value used in the merging is
#'   plotted at each node.
#' @param isCount logical as to whether input data is a count matrix. See
#'   details.
#' @param doPlot logical as to whether to plot the dendrogram (overrides 
#'   \code{plotType} value). Mainly used for internal coding purposes.
#' @param dendroSamples If x is a matrix, this is a dendrogram on the samples
#'   (unlike \code{dendro} which is a dendrogram on the clusters); this should
#'   be a dendrogram that is the same topology as the dendrogram in
#'   \code{dendro}, but includes individual entries for the samples (see
#'   \code{\link{makeDendrogram}}). This is used ONLY for plotting the
#'   clusterings before and after merging (if \code{plotType} is not 'none'). If
#'   x is a \code{ClusterExperiment} object, this is passed internally and is
#'   not specified by the user.
#' @param ... for signature \code{matrix}, arguments passed to the 
#'   \code{\link{plot.phylo}} function of \code{ade4} that plots the dendrogram.
#'   For signature \code{ClusterExperiment} arguments passed to the method for 
#'   signature \code{matrix} and then onto \code{\link{plot.phylo}}.
#' @inheritParams clusterMany,matrix-method
#'
#' @details If  \code{isCount=TRUE}, and the input is a matrix,
#'   \code{log2(count + 1)} will be used for \code{\link{makeDendrogram}} and the
#'   original data with voom correction will be used in
#'   \code{\link{getBestFeatures}}). If input is
#'   \code{\link{ClusterExperiment}}, then setting \code{isCount=TRUE} also means
#'   that the log2(1+count) will be used as the transformation, like for
#'   the matrix case as well as the voom calculation, and will NOT use the
#'   transformation stored in the object. If FALSE, then transform(x) will be
#'   given to the input and will be used for both \code{makeDendrogram} and
#'   \code{getBestFeatures}, with no voom correction.
#' @details "JC" refers to the method of Ji and Cai (2007), and implementation
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
#'   indicate where the clusters will be merged (assuming \code{plotType} is not 'none').
#' @return If `x` is a matrix, it returns (invisibly) a list with elements
#'   \itemize{ \item{\code{clustering}}{ a vector of length equal to ncol(x)
#'   giving the integer-valued cluster ids for each sample. "-1" indicates the
#'   sample was not clustered.} \item{\code{oldClToNew}}{ A table of the old
#'   cluster labels to the new cluster labels.} \item{\code{propDE}}{ A table of
#'   the proportions that are DE on each node.}
#'   \item{\code{originalClusterDendro}}{ The dendrogram on which the merging
#'   was based (based on the original clustering).} }
#' @return If `x` is a \code{\link{ClusterExperiment}}, it returns a new
#'   \code{ClusterExperiment} object with an additional clustering based on the
#'   merging. This becomes the new primary clustering.
#' @examples
#' data(simData)
#'
#' #create a clustering, for 8 clusters (truth was 3)
#' cl<-clusterSingle(simData, clusterFunction="pam", subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))
#'
#' #make dendrogram
#' cl <- makeDendrogram(cl)
#'
#' #merge clusters with plotting. Note argument 'use.edge.length' can improve
#' #readability
#' merged <- mergeClusters(cl, plotType="all",
#' mergeMethod="adjP", use.edge.length=FALSE)
#'
#' #compare merged to original
#' table(primaryCluster(cl), primaryCluster(merged))
#' @export
#' @importFrom phylobase labels descendants ancestors getNode
#' @importClassesFrom phylobase phylo4
#' @importFrom graphics plot
#' @importFrom ape plot.phylo phydataplot
#' @importFrom howmany howmany lowerbound
#' @importFrom locfdr locfdr
#' @rdname mergeClusters
setMethod(f = "mergeClusters",
          signature = signature(x = "matrix"),
          definition = function(x, cl, dendro=NULL,
                          mergeMethod=c("none", "adjP", "locfdr", "MB", "JC"),
                          plotType=c("none", "all", "mergeMethod","adjP", "locfdr", "MB", "JC"), 
                          cutoff=0.1, doPlot=TRUE,
                          isCount=TRUE, dendroSamples=NULL, ...) {
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
  plotType <- match.arg(plotType)
  if(mergeMethod=="none" & plotType=="none") stop("mergeMethod and plotType both equal 'none'; nothing to be done.")
  if(plotType=="mergeMethod" & mergeMethod=="none") {
    stop("can only plot merge method values if one method is selected")
  }

  #get test-statistics for the contrasts corresponding to each node (and return all)
  sigTable <- getBestFeatures(x, cl, contrastType=c("Dendro"), dendro=dendro,
                               contrastAdj=c("All"),
                              number=nrow(x), p.value=1, isCount=isCount)

  #divide table into each node.
  whMethodCalculate<-if(!mergeMethod=="none") mergeMethod else c()
  if(plotType=="all") whMethodCalculate<-c("adjP", "locfdr", "MB", "JC")
  if(plotType%in% c("adjP", "locfdr", "MB", "JC")) whMethodCalculate<-unique(c(whMethodCalculate,plotType))
  sigByNode <- by(sigTable, sigTable$ContrastName, function(x) {
      mb <-if("MB" %in% whMethodCalculate)  .myTryFunc(pvalues=x$P.Value, FUN=.m1_MB) else NA
      locfdr <-if("locfdr" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_locfdr) else NA
      jc <-if("JC" %in% whMethodCalculate)  .myTryFunc(tstats=x$t, FUN=.m1_JC) else NA
      adjP<-if("adjP" %in% whMethodCalculate)  .m1_adjP(x$adj) else NA
      return(c("adjP"=adjP, "locfdr"=locfdr, "MB"=mb,"JC"=jc))
  })
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

  #browser()
  nodePropTable <- do.call("rbind", sigByNode)
  annotTable <- data.frame("Node"=names(sigByNode),
                              "Contrast"=sigTable$Contrast[match(names(sigByNode), sigTable$ContrastName)])
    #add merge information:
   if(mergeMethod != "none" && length(whToMerge)>0 && length(which(whToMerge)) > 0){
     logicalMerge<-annotTable$Node %in%nodesToMerge

 } else logicalMerge<-rep(FALSE,length=nrow(annotTable))
  nodePropTable<-data.frame(annotTable,"Merged"=logicalMerge,nodePropTable)
  
  if(mergeMethod=="none"){
    newcl<-NULL #was just the original and nothing changed, so don't return something that makes it look like theres a new clustering
    oldClToNew<-NULL
      }
  else oldClToNew=table(Original=cl, New=newcl)
  out<-list(clustering=newcl, oldClToNew=oldClToNew,
                 propDE=nodePropTable, originalClusterDendro=dendro,mergeMethod=mergeMethod)
  if(doPlot){
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
	  if(is.null(dendroSamples)){
		  clMat<-unique(clMat)
		  rownames(clMat)<-as.character(clMat[,1])
	  }
	  #browser()
	  if(!is.null(dendroSamples)) .plotDendro(dendroSamples,leafType="samples",mergeOutput=out,mergePlotType=plotType,mergeMethod=mergeMethod,cl=clMat,label="name",outbranch=any(cl<0),...)
	 else .plotDendro(dendro,leafType="clusters",mergeOutput=out,mergePlotType=plotType,mergeMethod=mergeMethod,cl=clMat,label="name",...)
  	
  }
  invisible(out)
}
)


#' @rdname mergeClusters
#' @export
#' @param clusterLabel a string used to describe the type of clustering. By 
#'   default it is equal to "mergeClusters", to indicate that this clustering is
#'   the result of a call to mergeClusters.
#' @param labelLeaves if plotting, then whether leaves of dendrogram should be
#'   labeled by rectangular blocks of color ("colorblock")  or with the names of
#'   the leaves ("name").
#' @param leaves if plotting, whether the leaves should be the clusters or the
#'   samples. Choosing 'samples' allows for visualization of how many samples.
#' @details Note that \code{leaves='samples'} is currently fragile, in the sense
#'   that the alignment of the nodes in the cluster dendrogram (which correspond
#'   to the merge cutoff values) to that of the dendrogram with individual
#'   sample values is fragile, and may not be correct.
setMethod(f = "mergeClusters",
          signature = signature(x = "ClusterExperiment"),
          definition = function(x, eraseOld=FALSE,isCount=FALSE,
                                mergeMethod="none",plotType="all",clusterLabel="mergeClusters",leaves=c("clusters","samples" ),labelLeaves=c("name","colorblock","ids"),...) {
  labelLeaves<-match.arg(labelLeaves)
  leaves<-match.arg(leaves)
  if(is.null(x@dendro_clusters)) {
    stop("`makeDendrogram` needs to be called before `mergeClusters`")
  }
  else{
    cl<-clusterMatrix(x)[,x@dendro_index]
    note("Merging will be done on '",clusterLabels(x)[x@dendro_index],"', with clustering index",x@dendro_index)
  }
  if(isCount) note("If `isCount=TRUE` the data will be transformed with voom() rather than
with the transformation function in the slot `transformation`.
This makes sense only for counts.")
  
###Note, doPlot=FALSE, and then manually call .plotDendro afterwards to allow for passage of colors, etc.
  outlist <- mergeClusters(x=if(!isCount) transform(x) else assay(x),
                           cl=cl,
                           dendro=x@dendro_clusters, plotType=plotType,doPlot=FALSE,
                           isCount=isCount,mergeMethod=mergeMethod, ...)
  
  if(mergeMethod!="none"){#only add a new cluster if there was a mergeMethod. otherwise, mergeClusters just returns original cluster!
    #----
    #add "m" to name of cluster
    #----
    newObj <- clusterExperiment(x, outlist$clustering,
                                transformation=transformation(x),
                                clusterTypes="mergeClusters")
    #add "m" to name of cluster
    newObj<-.addPrefixToClusterNames(newObj,prefix="m",whCluster=1)
    clusterLabels(newObj) <- clusterLabel
    ##Check if pipeline already ran previously and if so increase
    x<-.updateCurrentWorkflow(x,eraseOld,"mergeClusters")
    if(!is.null(x)) retval<-.addNewResult(newObj=newObj,oldObj=x)
    else retval<-.addBackSEInfo(newObj=newObj,oldObj=x)
  }
  else{ #don't do anything, since there was no merging done.
    retval<-x
  }
  if(plotType!="none"){
      dend<- switch(leaves,"samples"=retval@dendro_samples,"clusters"=retval@dendro_clusters)
  		leg<-clusterLegend(retval)[[retval@dendro_index]]
      cl<-switch(leaves,"samples"=clusterMatrix(retval)[,retval@dendro_index],"clusters"=NULL)
  	if(leaves=="samples") names(cl)<-colnames(retval)
      if(labelLeaves=="id") leg[,"name"]<-leg[,"clusterIds"]
  	label<-switch(labelLeaves,"name"="name","colorblock"="colorblock","ids"="name")
  	outbranch<-FALSE
  	if(leaves=="samples" & any(cl<0)) outbranch<-TRUE

  # outbranch<-any(clusterMatrix(retval)[,retval@dendro_index]<0)
  # cl<-clusterMatrix(retval,whichCluster=retval@dendro_index)
  # rownames(cl)<-colnames(retval)
  # dend<-ifelse(leaves=="samples", retval@dendro_samples,retval@dendro_clusters)
     .plotDendro(dendro=dend,leafType=leaves,mergeOutput=outlist,mergePlotType=plotType,mergeMethod=mergeMethod,cl=cl,clusterLegendMat=leg,label=label,outbranch=outbranch)

 # .plotDendro(retval@dendro_clusters,leafType="clusters",mergeOutput=outlist,mergePlotType=plotType,mergeMethod=mergeMethod,cl=clusterMatrix(retval,whichCluster=retval@dendro_index),clusterLegendMat=clusterLegend(retval)[[retval@dendro_index]],label="colorblock")
  }
  
  invisible(retval)
}
)


.plotDendro<-function(dendro,leafType="clusters",mergePlotType=NULL,mergeMethod=NULL,mergeOutput=NULL,clusterLegendMat=NULL,cl=NULL,label=c("name","colorblock"),outbranch=FALSE,...){
	label<-match.arg(label)
    phylo4Obj <- .makePhylobaseTree(dendro, "dendro",isSamples=(leafType=="samples"),outbranch=outbranch)
    phyloObj <- as(phylo4Obj, "phylo")
	#browser()
	plotArgs<-list(...)
	###############
	### For plotting of dendrogram for the merging
	### Add information about the merging
	###############
	if(!is.null(mergePlotType) && mergePlotType %in% c("all","adjP", "locfdr", "MB", "JC","mergeMethod")){
        #####
        #convert names of internal nodes for plotting
        #####
        #match to order of tree
		#browser()
	    sigInfo<-mergeOutput$propDE
	    whToMerge<-which(sigInfo$Merged)
	    nodesToMerge<-sigInfo$Node[whToMerge]
	    methods<-colnames(sigInfo[,-c(1:3)])
        m <- match( sigInfo$Node,phyloObj$node)
		if(any(is.na(m))) stop("some nodes in mergeOutput not in the given dendrogram")
        edgeLty <- rep(1, nrow(phyloObj$edge))
        if(mergeMethod != "none" && length(whToMerge) > 0) {
            #which of nodes merged
			whMerge <- which(phyloObj$node.label %in% nodesToMerge) 
            nodeNumbers <- (length(phyloObj$tip) + 1):max(phyloObj$edge)
            whEdge <- which(phyloObj$edge[,1] %in% nodeNumbers[whMerge])
            edgeLty[whEdge] <- 2
        }
        if(mergePlotType == "mergeMethod"){
            if(!mergeMethod %in% methods) stop("mergeMethod not in methods of output")
            phyloObj$node.label[m] <- as.character(signif(sigInfo[,mergeMethod],2))
        }
        if(mergePlotType %in% c("all","adjP", "locfdr", "MB", "JC")) {
            meth<-if(mergePlotType=="all") methods else methods[methods%in%mergePlotType]
            phyloObj$node.label 
			phyloObj$node.label[m] <- apply(sigInfo[,meth,drop=FALSE],1, function(x){
                whKp<-which(!is.na(x))
                paste(paste(meth[whKp], signif(x[whKp],2), sep=":"), collapse=",\n")})
        }
		phyloObj$node.label[-m]<-""
		plotArgs$show.node.label<-TRUE
		plotArgs$edge.lty<-edgeLty
	}
	###############
	### Generic:
    ### Add color of cluster and cluster/sample name from the object.
	###############
	#temporary, do only 1 clustering:
	if(is.matrix(cl) && ncol(cl)>1) cl<-cl[,1,drop=FALSE]
	if(label=="colorblock" & is.null(clusterLegendMat)){
		#create a default color scheme
		clusterIds<-sort(unique(cl))
		clusterLegendMat<-cbind("clusterIds"=clusterIds,"name"=clusterIds,"color"=bigPalette[1:length(clusterIds)])
	}
    if(!is.null(clusterLegendMat)){
		if(leafType=="clusters"){
			m<-match(phyloObj$tip.label,clusterLegendMat[,"clusterIds"])
	        if(any(is.na(m))) stop("clusterIds do not match dendrogram labels")
	        phyloObj$tip.label<-clusterLegendMat[m,"name"]
	        tip.color<-clusterLegendMat[m,"color"]
			if(label=="colorblock"){
				#browser()
				clusterLegendMat<-clusterLegendMat[!clusterLegendMat[,"clusterIds"]%in%c(-1,-2),]
				colorMat<-matrix(clusterLegendMat[,"name"],ncol=1)
				row.names(colorMat)<-clusterLegendMat[,"name"]
				cols<-clusterLegendMat[,"color"]
				names(cols)<-clusterLegendMat[,"name"]
				
				#code that actually maps to the colors:
			    # lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
			    # x <- .matchDataPhylo(x, phy)
			    # n <- length(phy$tip.label)
				# one2n <- seq_len(n)
				# y1 <- lastPP$yy[one2n]
			    # o <- order(y1)
	            # x <- if (style == "image") x[o, o]
	            # else if (is.vector(x)) x[o]
	            # else x[o, ]
				#nux <- length(ux <- unique.default(x))
				#x <- match(x, ux)
				#co <- funcol(nux)
				#rect(xl, yb, xr, yt, col = co[x], xpd = TRUE, ...)
				# so colors need to be in the order of unique.default(x)
			}
		
		}
		else{
			m<-match(cl,clusterLegendMat[,"clusterIds"])
	        tip.color<-clusterLegendMat[m,"color"]		
			if(label=="colorblock"){
				colorMat<-matrix(clusterLegendMat[m,"name"],ncol=1)
				rownames(colorMat)<-names(cl)
				cols<-tip.color
				names(cols)<-clusterLegendMat[m,"name"]
				
			}	
		}
    }
    else tip.color<-"black"
		
	###############
	#this next code is hack to deal with error sometimes get if very long edge length -- usually due to unusual distance, etc.
	# Divides edge lengths so not too large.
	###############
	if(max(phyloObj$edge.length)>1e6) phyloObj$edge.length <- phyloObj$edge.length / max(phyloObj$edge.length) 
		
		
	
	#	browser()
	if(label=="name") do.call(ape::plot.phylo,c(list(phyloObj, tip.color=tip.color),plotArgs))
	else{#if colorblock
		phyloPlotOut<-do.call(ape::plot.phylo,c(list(phyloObj, tip.color=tip.color,show.tip.label=FALSE,plot=FALSE),plotArgs))
		treeWidth<-phyloPlotOut$x.lim[2]
		do.call(ape::plot.phylo,c(list(phyloObj, tip.color=tip.color,show.tip.label=FALSE,x.lim=treeWidth*1.5),plotArgs))
		#this is a temporary hack, because right now function has bug and fails for a 1-column matrix or vector. Have reported this 5/23/2017.
		if(ncol(colorMat)==1){
			colorMat<-cbind(colorMat,colorMat)
		}
		
		
		#we have to do this to get order for colors to be what we want!
		#basically have to redo code in phydataplot so figure out what order is in plot of the leaves, etc. Poor function. 
		#this doesn't work! can't find .PlotPhyloEnv 
		# added ape:::, perhaps will work. But don't know how I can export it in package???
		getColFun<-function(x,phy,namedColors){
			x <- ape:::.matchDataPhylo(x, phy)
			n <- length(phy$tip.label)
			one2n <- seq_len(n)
			lastPP <- get("last_plot.phylo", envir = ape:::.PlotPhyloEnv)
			y1 <- lastPP$yy[one2n]
			o <- order(y1)
			ux<-unique.default(x[o])
			m<-match(as.character(ux),names(namedColors))
			function(n){namedColors[m]}
		}
		#browser()
		ape::phydataplot(x=colorMat, phy=phyloObj, style="mosaic",offset=treeWidth*.5/16, width = treeWidth*.5/4, border = NA, lwd = 3,legend = "side", funcol = getColFun(colorMat,phyloObj,cols))

		
	}
	
	invisible(phyloObj)
}


.myTryFunc<-function(FUN,...){
  x<-try(FUN(...),silent=TRUE)
  if(!inherits(x, "try-error")) return(x)
  else return(NA)
}

#functions for estimating m1/m, the proportion of non-null
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
