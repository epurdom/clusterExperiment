.availMergeMethods<-c("adjP", "locfdr", "MB", "JC")	
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
#'   was based (based on the original clustering).} }
#' @return If `x` is a \code{\link{ClusterExperiment}}, it returns a new 
#'   \code{ClusterExperiment} object with an additional clustering based on the 
#'   merging. This becomes the new primary clustering.
#' @seealso makeDendrogram, plotDendrogram, getBestFeatures
#' @examples
#' data(simData)
#' 
#' #create a clustering, for 8 clusters (truth was 3)
#' cl<-clusterSingle(simData, clusterFunction="pam", subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))
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
                          mergeMethod=c("none", "adjP", "locfdr", "MB", "JC"),
                          plotInfo=c("none", "all", "mergeMethod","adjP", "locfdr", "MB", "JC"), 
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

  #get test-statistics for the contrasts corresponding to each node (and return all)
  sigTable <- getBestFeatures(x, cl, contrastType=c("Dendro"), dendro=dendro,
                               contrastAdj=c("All"),
                              number=nrow(x), p.value=1, isCount=isCount)
#browser()
  #divide table into each node.
  whMethodCalculate<-if(!mergeMethod=="none") mergeMethod else c()
  if(plotInfo=="all") whMethodCalculate<-.availMergeMethods
  if(plotInfo%in% .availMergeMethods) whMethodCalculate<-unique(c(whMethodCalculate,plotInfo))
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
	  #browser()
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
                                mergeMethod="none",plotInfo="all",clusterLabel="mergeClusters",leafType=c("samples","clusters" ),labelType=c("colorblock","name","ids"),plot=TRUE,...) {
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
  outlist <- mergeClusters(x=if(!isCount) transform(x) else assay(x),
                           cl=cl,
                           dendro=x@dendro_clusters, plotInfo=plotInfo,plot=FALSE,
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
  if(plot){
    dend<- switch(leafType,"samples"=retval@dendro_samples,"clusters"=retval@dendro_clusters)
  	# leg<-clusterLegend(retval)[[retval@dendro_index]]
  	#     cl<-switch(leafType,"samples"=clusterMatrix(retval)[,retval@dendro_index],"clusters"=NULL)
	#browser()
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
  #browser()
    if(labelType=="id") leg[,"name"]<-leg[,"clusterIds"]
  	label<-switch(labelType,"name"="name","colorblock"="colorblock","ids"="name")
  	outbranch<-FALSE
  	if(leafType=="samples" & any(cl<0)) outbranch<-retval@dendro_outbranch
		#if(leafType=="samples" & any(cl<0)) outbranch<-TRUE

  # outbranch<-any(clusterMatrix(retval)[,retval@dendro_index]<0)
  # cl<-clusterMatrix(retval,whichCluster=retval@dendro_index)
  # rownames(cl)<-colnames(retval)
  # dend<-ifelse(leafType=="samples", retval@dendro_samples,retval@dendro_clusters)
     .plotDendro(dendro=dend,leafType=leafType,mergeOutput=outlist,mergePlotType=plotInfo,mergeMethod=mergeMethod,cl=cl,clusterLegendMat=leg,label=label,outbranch=outbranch,removeOutbranch=outbranch)
  }
  
  invisible(retval)
}
)




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
