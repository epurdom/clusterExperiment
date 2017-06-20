#' @title Plot dendrogram of clusterExperiment object
#'   
#' @description Plots the dendrogram saved in a clusterExperiment object
#'   
#' @param x a \code{\link{ClusterExperiment}} object.
#' @param leafType if "samples" the dendrogram has one leaf per sample,
#'   otherwise it has one per cluster.
#' @param main passed to the \code{plot.phylo} function to set main title.
#' @param sub passed to the \code{plot.phylo} function to set subtitle.
#' @param labelType one of 'name', 'colorblock' or 'id'. If 'Name' then 
#'   dendrogram will be plotted, and name of cluster or sample (depending on 
#'   type of value for \code{leafType}) will be plotted next to the leaf of the 
#'   dendrogram. If 'colorblock', rectangular blocks, corresponding to the color
#'   of the cluster will be plotted, along with cluster name legend. If 'id' the
#'   internal clusterIds value will be plotted (only appropriate if 
#'   \code{leafType="clusters"}).
#' @param ... arguments passed to the \code{\link{plot.phylo}} function of
#'   \code{ape} that plots the dendrogram.
#' @param whichClusters only used if \code{leafType="samples"}). If numeric, an
#'   index for the clusterings to be plotted with dendrogram. Otherwise, 
#'   \code{whichClusters} can be a character value identifying the 
#'   \code{clusterTypes} to be used, or if not matching \code{clusterTypes} then
#'   \code{clusterLabels}; alternatively \code{whichClusters} can be either 
#'   'all' or 'workflow' or 'primaryCluster' to indicate choosing all clusters
#'   or choosing all \code{\link{workflowClusters}}. Default 'dendro' indicates
#'   using the clustering that created the dendrogram.
#'   @param removeOutgroup logical, only applicable if there are missing samples
#'     (i.e. equal to -1 or -2), \code{leafType="samples"} and the dendrogram
#'     for the samples was made by putting missing samples in an outbranch. In
#'     which case, if this parameter is TRUE, the outbranch will not be plotted,
#'     and if FALSE it will be plotted.
#' @param legend logical, only applicable if \code{labelType="colorblock"}.
#'   Passed to \code{\link{phydataplot}} in \code{\link{ape}} package that is
#'   used to draw the color values of the clusters/samples next to the
#'   dendrogram. Options are 'none', 'below', or 'side'
#' @aliases plotDendrogram
#' @details If \code{leafType="clusters"}, the plotting function will work best
#'   if the clusters in the dendrogram correspond to the primary cluster. This
#'   is because the function colors the cluster labels based on the colors of
#'   the clusterIds of the primaryCluster
#' @importFrom ape plot.phylo
#' @export
#' 
#' @examples
#' data(simData)
#' 
#' #create a clustering, for 8 clusters (truth was 3) 
#' cl <-clusterSingle(simData, subsample=FALSE, 
#' sequential=FALSE, clusterDArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#' 
#' #create dendrogram of clusters and then 
#' # merge clusters based ondendrogram: 
#' cl <- makeDendrogram(cl) 
#' cl <- mergeClusters(cl,mergeMethod="adjP",cutoff=0.1,plot=FALSE) 
#' plotDendrogram(cl) 
#' plotDendrogram(cl,leafType="samples",whichClusters="all",labelType="colorblock")
#' 
#' @export
#' @rdname plotDendrogram
setMethod(
  f = "plotDendrogram",
  signature = "ClusterExperiment",
  definition = function(x,whichClusters="dendro",leafType=c("clusters","samples" ),  labelType=c("name","colorblock","ids"), main,sub,removeOutbranch=TRUE,legend='side',...)
  {
    if(is.null(x@dendro_samples) || is.null(x@dendro_clusters)) stop("No dendrogram is found for this ClusterExperiment Object. Run makeDendrogram first.")
    leafType<-match.arg(leafType)
	labelType<-match.arg(labelType)
	whCl<-.TypeIntoIndices(x,whClusters=whichClusters)
    if(length(whCl)==0) stop("given whichClusters value does not match any clusters")

    if(missing(main)) main<-ifelse(leafType=="samples","Dendrogram of samples", "Dendrogram of clusters")
    if(missing(sub)) sub<-paste("Dendrogram made with '",clusterLabels(x)[dendroClusterIndex(x)],"', cluster index ",dendroClusterIndex(x),sep="")
	dend<- switch(leafType,"samples"=x@dendro_samples,"clusters"=x@dendro_clusters)

    cl<-switch(leafType,"samples"=clusterMatrix(x)[,whCl,drop=FALSE],"clusters"=NULL)
	if(leafType=="samples") rownames(cl)<-if(!is.null(colnames(x))) colnames(x) else as.character(1:ncol(x))
	if(length(whCl)==1){
		leg<-clusterLegend(x)[[whCl]]
	    if(labelType=="id") leg[,"name"]<-leg[,"clusterIds"]		
	}
	else{
		leg<-clusterLegend(x)[whCl]
	    if(labelType=="id") leg<-lapply(leg,function(x){x[,"name"]<-x[,"clusterIds"]; return(x)})	
	}
	label<-switch(labelType,"name"="name","colorblock"="colorblock","ids"="name")
	invisible(.plotDendro(dendro=dend,leafType=leafType,mergeMethod=NULL,mergeOutput=NULL,clusterLegendMat=leg,cl=cl,label=label,outbranch=x@dendro_outbranch,main=main,sub=sub,removeOutbranch=removeOutbranch,legend=legend,...))
    
  })
  



########
# Internal plotting function used by both mergeClusters and plotDendrogram
#' @importFrom phylobase labels descendants ancestors getNode edgeLength rootNode nodeLabels nNodes subset
#' @importClassesFrom phylobase phylo4 
#' @importFrom graphics plot
#' @importFrom ape plot.phylo phydataplot
.plotDendro<-function(dendro,leafType="clusters",mergePlotType=NULL,mergeMethod=NULL,mergeOutput=NULL,clusterLegendMat=NULL,cl=NULL,label=c("name","colorblock"),outbranch=FALSE,removeOutbranch=FALSE,legend="below",...){
	label<-match.arg(label)
	phylo4Obj <- .makePhylobaseTree(dendro, "dendro",isSamples=(leafType=="samples"),outbranch=outbranch)
	#---
	#remove the outbranch from the dendrogram and from cl
	#(note this is using phylo4 obj)
	#---
	if(outbranch & removeOutbranch & leafType=="samples"){
		rootNode<-phylobase::rootNode(phylo4Obj)
		rootChild<-phylobase::descendants(phylo4Obj,node=rootNode,type="children")
		tips<-phylobase::getNode(phylo4Obj,type="tip")
		whMissingNode<-grep("MissingNode",names(rootChild))
		if(length(whMissingNode)==0){
			#check not a single -1 sample from root:
			if(any(rootChild %in% tips)){
        #which ever rootChild is in tips must be single missing sample 
			  #because can't make dendrogram with only 1 cluster so couldn't run plot or mergeClusters. 
			  clusterNode<-rootChild[!rootChild %in% tips]
			  #stop("Internal coding error: need to fix .plotDendro to deal with when single missing sample")
			}
			else stop("Internal coding error: no outbranch nodes")	
		} 
		else clusterNode<-rootChild[-whMissingNode]
		if(length(clusterNode)!=1) stop("Internal coding error: removing missing node does not leave exactly 1 descendent of root")
		clusterTips<-phylobase::descendants(phylo4Obj,node=clusterNode,type="tip")
		if(length(clusterTips)==0) stop("Internal coding error: no none missing samples in tree")
		namesClusterTips<-names(clusterTips)
		if(is.matrix(cl)) cl<-cl[namesClusterTips,] else cl<-cl[namesClusterTips]
		phylo4Obj<-phylobase::subset(phylo4Obj, node.subtree=clusterNode)
		#set outbranch=FALSE because now doesn't exist in tree...
		outbranch<-FALSE
	}
	phyloObj <- as(phylo4Obj, "phylo")

  	plotArgs<-list(...)
	dataPct<-0.5
	offsetDivide<-16
	if(label=="colorblock" && is.null(cl) && leafType=="samples") stop("Internal coding error: must provide a clustering if label='colorblock'")
  	###############
  	### For plotting of dendrogram for the merging
  	### Add information about the merging as node labels and change edge type
  	###############
  	if(!is.null(mergePlotType) && mergePlotType %in% c("all","adjP", "locfdr", "MB", "JC","mergeMethod")){
          #####
          #convert names of internal nodes for plotting
          #####
          #match to order of tree
  		#browser()
  	    sigInfo<-mergeOutput$propDE
  	    whToMerge<-which(sigInfo$Merged)
  	    nodesToMerge<-as.character(sigInfo$Node[whToMerge])
  	    methods<-colnames(sigInfo[,-c(1:3)])
        m <- match( as.character(sigInfo$Node),phyloObj$node)
  		if(any(is.na(m))) stop("some nodes in mergeOutput not in the given dendrogram")
        edgeLty <- rep(1, nrow(phyloObj$edge))
		if(mergeMethod != "none" && length(whToMerge) > 0){
		  #which of nodes merged
		  whMerge <- which(phyloObj$node.label %in% nodesToMerge) 
		  nodeNumbers <- (length(phyloObj$tip) + 1):max(phyloObj$edge)
		  whEdge <- which(phyloObj$edge[,1] %in% nodeNumbers[whMerge])
		  edgeLty[whEdge] <- 2
		}
		if(mergePlotType == "mergeMethod"){
		  if(!mergeMethod %in% methods) stop("mergeMethod not in methods of output")
		  valsNodes<-as.character(signif(sigInfo[,mergeMethod],2))
		  valsNodes[is.na(valsNodes)]<-"NA" #make them print out as NA -- otherwise doesn't plot
		  phyloObj$node.label[m] <- valsNodes
		  # offsetDivide<-3
		  # dataPct<-.7
		}
		if(mergePlotType %in% c("all",.availMergeMethods)) {
		  	meth<-if(mergePlotType=="all") methods else methods[methods%in%mergePlotType]
		  	phyloObj$node.label[m] <- apply(sigInfo[,meth,drop=FALSE],1, 
			  function(x){
		  		whKp<-which(!is.na(x))
		  	  	vals<-paste(paste(meth[whKp], signif(x[whKp],2), sep=":"), collapse="\n")
				vals[is.na(vals)]<-"NA"
				return(vals)
				})
			if(mergePlotType!="all"){
			  # offsetDivide<-3
			  # dataPct<-.7
			}
			else{
			  # offsetDivide<-2.5
			  # dataPct<-.7

			}
		}
		  
  		phyloObj$node.label[-m]<-""
  		plotArgs$show.node.label<-TRUE
  		plotArgs$edge.lty<-edgeLty
  	}
  	###############
  	### Deal with clusterLegend object: 
	### - Make default if not provided and 
	### - If # of clusterings>1 make clusterLegend and cl matrix appropriate
  	###############
	if(label=="colorblock"){
		clusterLegend<-TRUE #doesn't do anything right now because phydataplot doesn't have option of no legend...
		if(is.null(clusterLegendMat)){ 
			#----
			#make default colors, works for vector or matrix cl
			#----
  			clusterIds<-sort(unique(as.vector(cl)))
			clusterLegendMat<-cbind("clusterIds"=clusterIds,"name"=clusterIds,"color"=bigPalette[1:length(clusterIds)])
  		}
		else{
			if(is.matrix(cl) && ncol(cl)>1){
			  	#if not provide list of cluster legends, do only 1st clustering provided (temporary while fixing so works for matrix)
				if(!is.list(clusterLegendMat) ) cl<-cl[,1,drop=FALSE]
				else{
					#----
					#create one big cl/clusterLegendMat object that will allow for coloring that is okay.
					#----
					nclusters<-ncol(cl)
					if(length(clusterLegendMat)!=nclusters) stop("Internal coding error -- wrong length of colors for clustering")
					newClusterLegendMat<-clusterLegendMat[[1]]
					newCl<-cl[,1]
					#make it general in case some day want more than just 2 clusterings
					for(ii in 2:nclusters){
						#browser()
						currMat<-clusterLegendMat[[ii]]
						currCl<-cl[,ii]
						whExistingColor<-which(currMat[,"color"] %in% newClusterLegendMat[,"color"])
						
						if(length(whExistingColor)>0){
							#find new id to give it
							matchNew<-match(currMat[whExistingColor,"color"],newClusterLegendMat[,"color"])
							oldId<-currMat[whExistingColor,"clusterIds"]
							newId<-newClusterLegendMat[matchNew,"clusterIds"]
							mexist<-match(currCl,oldId)
							newFullId<-as.numeric(newId[mexist])
							currCl[!is.na(mexist)]<-newFullId[!is.na(mexist)]
								
							#change name so combination, if not already the same
							whDiff<-which(newClusterLegendMat[matchNew,"name"]!=currMat[whExistingColor,"name"])
							if(length(whDiff)>0){
								combName<-paste(newClusterLegendMat[matchNew,"name"],currMat[whExistingColor,"name"],sep="/")
								newClusterLegendMat[matchNew[whDiff],"name"]<-combName[whDiff]
								
							}
							
							#remove from current color scheme
							currMat<-currMat[-whExistingColor,,drop=FALSE]
						}
						#browser()
						if(nrow(currMat)>0){
							## increase remaing ids 
							maxNew<-max(as.vector(newCl))
							oldId2<-currMat[,"clusterIds"]
							newId2<-seq(from=maxNew+1,by=1,length=length(oldId2)) #replace with this in legend
							mexist2<-match(currCl,oldId2) #match old ids to the clusterings vector
							newFullId2<-as.numeric(newId2[mexist2]) #will get NAs for those that don't match (e.g. have been changed by previous step)
							currCl[!is.na(mexist2)]<-newFullId2[!is.na(mexist2)]
							
							## change ids in currMat
							currMat[,"clusterIds"]<-newId2
							#browser()
							
							## test correct that no overlap in ids or names or colors:
							if(any(currMat[,"clusterIds"] %in% newClusterLegendMat[,"clusterIds"])) stop("Internal coding error: still overlap in cluster Ids")
							if(any(currMat[,"color"] %in% newClusterLegendMat[,"color"])) stop("Internal coding error: still overlap in color")
							
							## add to new cluster color legend
							newClusterLegendMat<-rbind(newClusterLegendMat,currMat)
						}
						newCl<-cbind(newCl,currCl)
						
					}
					#browser()
					clusterLegendMat<-newClusterLegendMat
					colnames(newCl)<-colnames(cl)
					rownames(newCl)<-rownames(cl)
					cl<-newCl
					clusterLegend<-FALSE
					
				}
					
			}
		}
	} 
	###############
	### Deal with clusterLegend object: 
	### - Add color of cluster and cluster/sample name to tip labels if labelType=="name"
	### - Make colorMat matrix if labelType=="colorblock"
	###############
	edge.width=1
	if(!is.null(clusterLegendMat)){
		if(leafType=="clusters"){
			#get rid of matching string 
			m<-match(gsub("ClusterId","",phyloObj$tip.label),clusterLegendMat[,"clusterIds"])
			#browser()
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
	

			}

		}
		if(leafType=="samples"){
			if(is.matrix(cl) && ncol(cl)>1){
				clNames<-row.names(cl)
				if(label=="colorblock"){
					colorMat<-apply(cl,2,function(x){
						m<-match(x,clusterLegendMat[,"clusterIds"])
						clusterLegendMat[m,"name"]
					})
					if(any(dim(colorMat)!=dim(cl))) stop("Internal coding error: dimensions of colorMat don't match input")
					dimnames(colorMat)<-dimnames(cl)
					#m<-match(cl[,1],clusterLegendMat[,"clusterIds"])
			    	cols<-clusterLegendMat[,"color"]
					names(cols)<-clusterLegendMat[,"name"]
	
				}
				tip.color<-"black"
				#browser()
			}
			else{
				if(is.matrix(cl)) cl<-cl[,1]
				clNames<-names(cl)
				m<-match(cl,clusterLegendMat[,"clusterIds"])
			    tip.color<-clusterLegendMat[m,"color"]		
				if(label=="colorblock"){
					colorMat<-matrix(clusterLegendMat[m,"name"],ncol=1)
					rownames(colorMat)<-names(cl)
					cols<-tip.color
					names(cols)<-clusterLegendMat[m,"name"]
				}	
			}
			if(label=="colorblock"){
				ntips<-length(phyloObj$tip.label)
				whClusterNode<-which(!is.na(phyloObj$node.label) & phyloObj$node.label!="")+ ntips
				#only edges going to/from these nodes
				whEdgePlot<-which(apply(phyloObj$edge,1,function(x){any(x %in% whClusterNode)}))
				edge.width<-rep(0,nrow(phyloObj$edge))
				edge.width[whEdgePlot]<-1
			}
			m<-match(phyloObj$tip.label,clNames)
		    if(any(is.na(m))) stop("names of cl do not match dendrogram labels")
			
		}
	}
	else tip.color<-"black"
  	

		#browser()
  	###############
  	#this next code is hack to deal with error sometimes get if very long edge length -- usually due to unusual distance, etc.
  	# Divides edge lengths so not too large.
  	###############
  	if(max(phyloObj$edge.length)>1e6) phyloObj$edge.length <- phyloObj$edge.length / max(phyloObj$edge.length) 
		
	prohibitOptions<-c("tip.color","node.pos","edge.width")
	if(any(prohibitOptions %in% names(plotArgs))) stop("User cannot set following options to plot.phylo:",paste(prohibitOptions, collapse=","))
	plotArgs<-c(plotArgs,list(tip.color=tip.color,node.pos=2,edge.width=edge.width))	
  	#	browser()
  	if(label=="name") do.call(ape::plot.phylo,c(list(phyloObj),plotArgs))
  	else{#if colorblock
  		phyloPlotOut<-do.call(ape::plot.phylo,c(list(phyloObj,show.tip.label = FALSE,plot=FALSE),plotArgs))
  		treeWidth<-phyloPlotOut$x.lim[2]
  		do.call(ape::plot.phylo,c(list(phyloObj,show.tip.label = FALSE,x.lim=treeWidth*(1+dataPct)),plotArgs))
  		
  		nclusters<-ncol(colorMat)
		colnames(colorMat)<-NULL		
		if(nclusters==1 & packageVersion("ape")<'4.1.0.6'){
			#this is a temporary hack, because right now function has bug and fails for a 1-column matrix or vector. Have reported this 5/23/2017 and now fixed in new version of ape.
		  			colorMat<-cbind(colorMat,colorMat)
		  		}
				
  		#we have to do this to get order for colors to be what we want!
  		#basically have to redo code in phydataplot so figure out what order is in plot of the leaves, etc. Poor function. New version of ape fixes this.
  		getColFun<-function(x,phy,namedColors){
  			x <- ape:::.matchDataPhylo(x, phy)
  			n <- length(phy$tip.label)
  			one2n <- seq_len(n)
  			lastPP <- get("last_plot.phylo", envir = ape:::.PlotPhyloEnv)
  			y1 <- lastPP$yy[one2n]
  			o <- order(y1)
			if(!is.null(ncol(x))) ux<-unique.default(x[o,])
  			else ux<-unique.default(x[o])
  			m<-match(as.character(ux),names(namedColors))
			namedColors[m]
  		}
  		if(packageVersion("ape")<'4.1.0.6') cols<-getColFun(colorMat,phyloObj,cols)
		colInput<-function(n){cols}
		ape::phydataplot(x=colorMat, phy=phyloObj, style="mosaic",offset=treeWidth*dataPct/offsetDivide, width = treeWidth*dataPct/4, border = NA, lwd = 3,legend = legend, funcol = colInput)
		if(nclusters>1 & !is.null(colnames(cl))){
			xloc<-treeWidth+treeWidth*dataPct/offsetDivide+seq(from=0,by=treeWidth*dataPct/4,length=ncol(cl))
			ypos<-par("usr")[4]+0*diff(par("usr")[3:4])
			text(x=xloc,y=ypos,labels=colnames(cl),srt=45,xpd=NA,adj=c(0,0))
			
		}
		#browser()
		
  	}
	
  	invisible(phyloObj)
  }
