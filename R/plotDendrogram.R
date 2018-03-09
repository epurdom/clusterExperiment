#' @title Plot dendrogram of ClusterExperiment object
#'   
#' @description Plots the dendrogram saved in a ClusterExperiment object
#'   
#' @param x a \code{\link{ClusterExperiment}} object.
#' @param leafType if "samples" the dendrogram has one leaf per sample,
#'   otherwise it has one per cluster.
#' @param main passed to the \code{plot.phylo} function to set main title.
#' @param sub passed to the \code{plot.phylo} function to set subtitle.
#' @param plotType one of 'name', 'colorblock' or 'id'. If 'Name' then 
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
#' @param removeOutbranch logical, only applicable if there are missing samples
#'     (i.e. equal to -1 or -2), \code{leafType="samples"} and the dendrogram
#'     for the samples was made by putting missing samples in an outbranch. In
#'     which case, if this parameter is TRUE, the outbranch will not be plotted,
#'     and if FALSE it will be plotted.
#' @param legend character, only applicable if \code{plotType="colorblock"}.
#'   Passed to \code{\link{phydataplot}} in \code{\link{ape}} package that is
#'   used to draw the color values of the clusters/samples next to the
#'   dendrogram. Options are 'none', 'below', or 'side'. (Note 'none' is only 
#'   available for 'ape' package >= 4.1-0.6).
#' @param nodeColors named vector of colors to be plotted on a node in the dendrogram. 
#' Names should match the name of the node (calls \code{\link[ape]{nodelabels}}). 
#' @param clusterLabelAngle angle at which label of cluster will be drawn. Only
#'	 applicable if \code{plotType="colorblock"}.
#' @param mergeInfo What kind of information about merge to plot on dendrogram. If not 
#' equal to "none", will replicate the kind of plot that \code{\link{mergeClusters}} 
#' creates, and the input to \code{mergeInfo} corresponds to that of \code{plotInfo} in \code{mergeClusters}.
#' @return A dendrogram is plotted. Returns (invisibly) a list with elements
#' \itemize{
#' \item{\code{plottedObject}}{ the \code{phylo} object that is plotted.}
#' \item{\code{originalObject}}{ the \code{phylo} object before adjusting the node/tip labels. }
#' }
#' @aliases plotDendrogram
#' @details If \code{leafType="clusters"}, the plotting function will work best
#'   if the clusters in the dendrogram correspond to the primary cluster. This
#'   is because the function colors the cluster labels based on the colors of
#'   the clusterIds of the primaryCluster
#' @importFrom ape plot.phylo nodelabels
#' @seealso \code{\link{mergeClusters}},\code{\link[ape]{plot.phylo}},\code{\link[ape]{nodelabels}},\code{\link[ape]{tiplabels}}
#' @export
#' 
#' @examples
#' data(simData)
#' 
#' #create a clustering, for 8 clusters (truth was 3) 
#' cl <-clusterSingle(simData, subsample=FALSE, 
#' sequential=FALSE, mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#' 
#' #create dendrogram of clusters and then 
#' # merge clusters based ondendrogram: 
#' cl <- makeDendrogram(cl) 
#' cl <- mergeClusters(cl,mergeMethod="adjP",cutoff=0.1,plot=FALSE) 
#' plotDendrogram(cl) 
#' plotDendrogram(cl,leafType="samples",whichClusters="all",plotType="colorblock")
#' 
#' @export
#' @rdname plotDendrogram
setMethod(
  f = "plotDendrogram",
  signature = "ClusterExperiment",
  definition = function(x,whichClusters="dendro",leafType=c("samples","clusters" ),  plotType=c("colorblock","name","ids"), mergeInfo="none", main, sub, clusterLabelAngle=45, removeOutbranch=TRUE, legend=c("side","below", "none"),nodeColors=NULL,...)
  {
    if(is.null(x@dendro_samples) || is.null(x@dendro_clusters)) stop("No dendrogram is found for this ClusterExperiment Object. Run makeDendrogram first.")
    leafType<-match.arg(leafType)
    plotType<-match.arg(plotType)
	
	possibleMergeValues<-c("none", "all","mergeMethod",.availMergeMethods)
    if(!is.null(x@merge_nodeProp)){
		otherVals<-colnames(x@merge_nodeProp)[!colnames(x@merge_nodeProp)%in%c("Node","Contrast")]
		possibleMergeValues<-unique(possibleMergeValues,otherVals)
    	
    }
	mergeInfo<-match.arg(mergeInfo,possibleMergeValues)
	
	mergeInfo<-match.arg(mergeInfo)
	legend<-match.arg(legend)
    whCl<-.TypeIntoIndices(x,whClusters=whichClusters)
    if(length(whCl)==0) stop("given whichClusters value does not match any clusters")
    if(leafType=="clusters" && whCl!=x@dendro_index){
      warning("if leafType=='clusters', 'whichClusters' must match the clusters that created the dendrogram. Changing 'leafType' to 'samples'")
      leafType<-"samples"
    }
    if(leafType=="clusters" & length(whCl)>1) stop("If leafType equal to 'clusters' 'whichClusters' must be of length 1 (i.e. single cluster).")
    if(missing(main)) main<-ifelse(leafType=="samples","Dendrogram of samples", "Dendrogram of clusters")
    if(missing(sub)) sub<-paste("Dendrogram made with '",clusterLabels(x)[dendroClusterIndex(x)],"', cluster index ",dendroClusterIndex(x),sep="")
    dend<- switch(leafType,"samples"=x@dendro_samples,"clusters"=x@dendro_clusters)
    
    cl<-switch(leafType,"samples"=clusterMatrix(x)[,whCl,drop=FALSE],"clusters"=NULL)
    if(leafType=="samples") rownames(cl)<-if(!is.null(colnames(x))) colnames(x) else as.character(1:ncol(x))
    if(length(whCl)==1){
      leg<-clusterLegend(x)[[whCl]]
      if(plotType=="id") leg[,"name"]<-leg[,"clusterIds"]		
    }
    else{
      leg<-clusterLegend(x)[whCl]
      if(plotType=="id") leg<-lapply(leg,function(x){x[,"name"]<-x[,"clusterIds"]; return(x)})	
    }
    label<-switch(plotType,"name"="name","colorblock"="colorblock","ids"="name")
#   mergePlotType=NULL,mergeMethod=NULL,mergeOutput=NULL, 

	if(is.na(x@merge_dendrocluster_index)) mergeInfo<-"none"
	if(mergeInfo=="none"){
		mergeInfo<-NULL
		mergeMethod<-NULL
	}
	else if(!is.na(x@merge_dendrocluster_index)){ #has node prop info, even if no merge cluster
		if(mergeInfo=="mergeMethod"){
			if(is.na(x@merge_index)){
				warning("Cannot plot merge method because there is none. Plotting all")
				mergeInfo<-"all"
			}
			else{
				mergeMethod<-x@merge_method
			}
		}
		else if(mergeInfo==x@merge_method ) mergeMethod<-x@merge_method #only do dotted lines if matches saved merged method -- is this good default or no? not sure .plotDendro works otherwise from ClusterExperiment object
		else mergeMethod<-"none"
	} else{
		warning("There is no information about merging -- will ignore input to 'mergeInfo'")
	}

	phyloOut<-.plotDendro(dendro=dend,leafType=leafType,mergeMethod=mergeMethod,mergePlotType=mergeInfo,mergeOutput=nodeMergeInfo(x),clusterLegendMat=leg,cl=cl,plotType=label,outbranch=x@dendro_outbranch,main=main,sub=sub,removeOutbranch=removeOutbranch,legend=legend,clusterLabelAngle=clusterLabelAngle,...)
	
	if(!is.null(nodeColors)){
		if(is.null(names(nodeColors))) warning("Must give names to node colors, ignoring argument nodeColors")
		m<-match(names(nodeColors),phyloOut$originalObject$node.label)
		nodeColors<-nodeColors[!is.na(m)]
		m<-m[!is.na(m)]
		if(length(m)>0){
			ape::nodelabels(rep("",length(nodeColors)),m+length(phyloOut$originalObject$tip.label), frame = "c", bg = nodeColors,cex=1/par("cex"))
		}
		else{warning("No names of node colors match node name, ignoring argument nodeColors")}
	}
	invisible(phyloOut)
    
  })




########
# Internal plotting function used by both mergeClusters and plotDendrogram
#' @importFrom phylobase labels descendants ancestors getNode edgeLength rootNode nodeLabels nNodes subset
#' @importClassesFrom phylobase phylo4 
#' @importFrom graphics plot
#' @importFrom ape plot.phylo phydataplot
.plotDendro<-function(dendro,leafType="clusters",mergePlotType=NULL,mergeMethod=NULL,mergeOutput=NULL,clusterLegendMat=NULL,cl=NULL,plotType=c("name","colorblock"),outbranch=FALSE,removeOutbranch=FALSE,legend="below",clusterLabelAngle=45,...){
	plotType<-match.arg(plotType)
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
		#
	
		if(is.matrix(cl)) cl<-cl[namesClusterTips,] else cl<-cl[namesClusterTips]
		phylo4Obj<-phylobase::subset(phylo4Obj, node.subtree=clusterNode)
		#set outbranch=FALSE because now doesn't exist in tree...
		outbranch<-FALSE
	}
	phyloObj <- as(phylo4Obj, "phylo")
  	plotArgs<-list(...)
	dataPct<-0.5
	offsetDivide<-16
	if(plotType=="colorblock" && is.null(cl) && leafType=="samples") stop("Internal coding error: must provide a clustering if plotType='colorblock'")
	origPhylo<-phyloObj #so can return
	
  	###############
  	### For plotting of dendrogram for the merging
  	### Add information about the merging as node labels and change edge type
	### Note: could probably have used nodelabels function and avoided some of this
  	###############
	if(!is.null(mergeOutput)){
		annotNames<-c("Node","Contrast","isMerged", "mergeClusterId")
		methods<-colnames(mergeOutput)[!colnames(mergeOutput)%in%annotNames] #possible for which have proportion saved
	}
	if(!is.null(mergePlotType) && !is.null(mergeOutput) && mergePlotType %in% c("all",methods,"mergeMethod")){
          #####
          #convert names of internal nodes for plotting
          #####
          #match to order of tree
  	    whToMerge<-which(mergeOutput$isMerged)
  	    nodesToMerge<-as.character(mergeOutput$Node[whToMerge])
		m <- match( as.character(mergeOutput$Node),phyloObj$node)
  		if(any(is.na(m))) stop("some nodes in merge node info not in the given dendrogram")
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
		  valsNodes<-as.character(signif(mergeOutput[,mergeMethod],2))
		  valsNodes[is.na(valsNodes)]<-"NA" #make them print out as NA -- otherwise doesn't plot
		  phyloObj$node.label[m] <- valsNodes
		  # offsetDivide<-3
		  # dataPct<-.7
		}
		else if(all(mergePlotType %in% c("all",methods))) {
			meth<-if(mergePlotType=="all") methods else methods[methods%in%mergePlotType]
		  	phyloObj$node.label[m] <- apply(mergeOutput[,meth,drop=FALSE],1, 
			  function(x){
		  		whKp<-which(!is.na(x))
				#fix up name of FC methods:
				mm<-meth[whKp]
				mm<-sapply(strsplit(mm,"_"),function(u){if(length(u)==2) return(sprintf("%s(>%s)",u[1],u[2])) else return(u)})
		  	  	vals<-paste(paste(mm, signif(x[whKp],2), sep=":"), collapse="\n")
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
		###Only set these if user doesn't...
		if(!"show.node.label" %in% names(plotArgs)) plotArgs$show.node.label<-TRUE
  		if(!"edge.lty" %in% names(plotArgs)) plotArgs$edge.lty<-edgeLty
  	}
  	###############
  	### Deal with clusterLegend object: 
	### - Make default if not provided and 
	### - If # of clusterings>1 make clusterLegend and cl matrix appropriate
  	###############
	
	if(plotType=="colorblock"){
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
						currMat<-clusterLegendMat[[ii]]
						currCl<-cl[,ii]
						
						#note that because subset to those samples that are not -1/-2 on
						# cl[,1], may have entire clusters in other columns of cl that disappear but still in color matrix with no entry in cl
						#reduce down the currMat to accomodate that
						whExist<-which(as.numeric(currMat[,"clusterIds"]) %in% currCl)
						currMat<-currMat[whExist, ,drop=FALSE]
						
						whExistingColor<-which(currMat[,"color"] %in% newClusterLegendMat[,"color"])
						
						if(length(whExistingColor)>0){
							#-----------
							#reassign the cluster id to the one matching existing color id.
							# only should affect names, not color of newClusterLegendMat
							#-----------
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
						
						if(nrow(currMat)>0){
							#-----------
							# for remaining, add to newClusterLegendMat 
							# with ids increased to be unique
							#-----------
							maxNew<-max(as.vector(newCl))
							oldId2<-currMat[,"clusterIds"]
							newId2<-seq(from=maxNew+1,by=1,length=length(oldId2)) #replace with this in legend
							mexist2<-match(currCl,oldId2) #match old ids to the clusterings vector
							newFullId2<-as.numeric(newId2[mexist2]) #will get NAs for those that don't match (e.g. have been changed by previous step)
							currCl[!is.na(mexist2)]<-newFullId2[!is.na(mexist2)]
							
							## change ids in currMat
							currMat[,"clusterIds"]<-newId2
							## test correct that no overlap in ids or names or colors:
							if(any(currMat[,"clusterIds"] %in% newClusterLegendMat[,"clusterIds"])){
								stop("Internal coding error: still overlap in cluster Ids")
							}
							if(any(currMat[,"color"] %in% newClusterLegendMat[,"color"])){
								stop("Internal coding error: still overlap in color")
							} 
							
							## add to new cluster color legend
							newClusterLegendMat<-rbind(newClusterLegendMat,currMat)
						}
						newCl<-cbind(newCl,currCl)
						
					}
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
	### Deal with clusterLegend object Part II: 
	### - Add color of cluster and cluster/sample name to tip labels if plotType=="name"
	### - Make colorMat matrix if plotType=="colorblock"
	###############
	edge.width=1
	if(!is.null(clusterLegendMat)){
		if(leafType=="clusters"){
			#get rid of matching string 
			m<-match(gsub("ClusterId","",phyloObj$tip.label),clusterLegendMat[,"clusterIds"])
		    if(any(is.na(m))) stop("clusterIds in clusterLegend do not match dendrogram labels")
		    phyloObj$tip.label<-clusterLegendMat[m,"name"]
		    tip.color<-clusterLegendMat[m,"color"]
			if(plotType=="colorblock"){
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
				if(plotType=="colorblock"){
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
			}
			else{
				if(is.matrix(cl)) cl<-cl[,1]
				clNames<-names(cl)
				m<-match(cl,clusterLegendMat[,"clusterIds"])
			    tip.color<-clusterLegendMat[m,"color"]		
				if(plotType=="colorblock"){
					colorMat<-matrix(clusterLegendMat[m,"name"],ncol=1)
					rownames(colorMat)<-names(cl)
					cols<-tip.color
					names(cols)<-clusterLegendMat[m,"name"]
				}	
			}
			if(plotType=="colorblock"){
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
	

  	#---------------
  	#this next code is hack to deal with error sometimes get if very long edge length -- usually due to unusual distance, etc.
  	# Divides edge lengths so not too large.
  	#---------------
  	if(max(phyloObj$edge.length)>1e6) phyloObj$edge.length <- phyloObj$edge.length / max(phyloObj$edge.length) 
		
	prohibitOptions<-c("tip.color","node.pos","edge.width","horizontal","direction","type")
	if(any(prohibitOptions %in% names(plotArgs))) stop("User cannot set following options to plot.phylo:",paste(prohibitOptions, collapse=","))
	plotArgs<-c(plotArgs,list(tip.color=tip.color,node.pos=2,edge.width=edge.width))	
  	if(plotType=="name") do.call(ape::plot.phylo,c(list(phyloObj),plotArgs))
  	else{
		#if colorblock
		#just calculate:
  	    if(!"show.tip.label"%in% plotArgs) plotArgs$show.tip.label<-FALSE
  		phyloPlotOut<-do.call(.calculatePlotPhylo,c(list(phyloObj,plot=FALSE),plotArgs))
  		treeWidth<-phyloPlotOut$x.lim[2]
		#plot dendrogram:
  		do.call(ape::plot.phylo,c(list(phyloObj,x.lim=treeWidth*(1+dataPct)),plotArgs))
  		
  		nclusters<-ncol(colorMat)
		colnames(colorMat)<-NULL		

		#########
		## Now require ape >5.0 so don't need this code. Not yet deleting...
		#########
		# if(nclusters==1 & packageVersion("ape")<'4.1.0.6'){
		# 	#this is a temporary hack, because right now function has bug and fails for a 1-column matrix or vector. Have reported this 5/23/2017 and now fixed in new version of ape.
		#   			colorMat<-cbind(colorMat,colorMat)
		#   		}				
		#   		#we have to do this to get order for colors to be what we want!
		#   		#basically have to redo code in phydataplot so figure out what order is in plot of the leaves, etc. Poor function. New version of ape fixes this.
		#   		getColFun<-function(x,phy,namedColors){
		#   			x <- ape:::.matchDataPhylo(x, phy)
		#   			n <- length(phy$tip.label)
		#   			one2n <- seq_len(n)
		#   			lastPP <- get("last_plot.phylo", envir = ape:::.PlotPhyloEnv)
		#   			y1 <- lastPP$yy[one2n]
		#   			o <- order(y1)
		# 	if(!is.null(ncol(x))) ux<-unique.default(x[o,])
		#   			else ux<-unique.default(x[o])
		#   			m<-match(as.character(ux),names(namedColors))
		# 	namedColors[m]
		#   		}
		# if(packageVersion("ape")<'4.1.0.6') cols<-getColFun(colorMat,phyloObj,cols)
		# if(packageVersion("ape")<'4.9.0.0' & packageVersion("ape")>='4.1.0.6') warning("If you have an ape package version between 4.1.0.6 and 5.0 installed, then the plotting of the clusters next to the dendrogram may not correctly use the colors in clusterLegend(object) nor have accurate legends." )
		colInput<-function(n){cols}
		width<-treeWidth*dataPct/nclusters
		ape::phydataplot(x=colorMat, phy=phyloObj, style="mosaic",offset=treeWidth*dataPct/offsetDivide, width = width, border = NA, lwd = 3,legend = legend, funcol = colInput)
		
		if(nclusters>1 & !is.null(colnames(cl))){
			xloc<-treeWidth+treeWidth*dataPct/offsetDivide+seq(from=0,by=width,length=nclusters)
			xloc<-xloc+width/2
			ypos<-par("usr")[4]-0.025*diff(par("usr")[3:4])	
			adj<-c(0,0)		
			if("cex" %in% names(list(...))) labcex<-list(...)[["cex"]]
				else labcex<-1	
			text(x=xloc,y=ypos,labels=colnames(cl),srt=clusterLabelAngle,xpd=NA,adj=adj,cex=labcex)
			
		}
		
  	}
	
  	invisible(list(plottedObject=phyloObj,originalObject=origPhylo))
  }


#' @importFrom ape node_height_clado node_height node_depth node_depth_edgelength is.ultrametric unrooted.xy
###Copy code from ape:plot.phylo so get how calculates x.lim. Just removed the plot parts.
.calculatePlotPhylo<-function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
      show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
      edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
      adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
      label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
      direction = "rightwards", lab4ut = NULL, tip.color = "black", 
      plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
      align.tip.label = FALSE, ...) 
  {
      Ntip <- length(x$tip.label)
      if (Ntip < 2) {
          warning("found less than 2 tips in the tree")
          return(NULL)
      }
      .nodeHeight <- function(edge, Nedge, yy) .C(ape::node_height, 
          as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
          as.double(yy))[[4]]
      .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(ape::node_depth, 
          as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
              2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
      .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
          edge.length) .C(ape::node_depth_edgelength, as.integer(edge[, 
          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
          double(Ntip + Nnode))[[5]]
      Nedge <- dim(x$edge)[1]
      Nnode <- x$Nnode
      if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
          stop("tree badly conformed; cannot plot. Check the edge matrix.")
      ROOT <- Ntip + 1
      type <- match.arg(type, c("phylogram", "cladogram", "fan", 
          "unrooted", "radial"))
      direction <- match.arg(direction, c("rightwards", "leftwards", 
          "upwards", "downwards"))
      if (is.null(x$edge.length)) {
          use.edge.length <- FALSE
      }
      else {
          if (use.edge.length && type != "radial") {
              tmp <- sum(is.na(x$edge.length))
              if (tmp) {
                  warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
                  use.edge.length <- FALSE
              }
          }
      }
      if (is.numeric(align.tip.label)) {
          align.tip.label.lty <- align.tip.label
          align.tip.label <- TRUE
      }
      else {
          if (align.tip.label) 
              align.tip.label.lty <- 3
      }
      if (align.tip.label) {
          if (type %in% c("unrooted", "radial") || !use.edge.length || 
              is.ultrametric(x)) 
              align.tip.label <- FALSE
      }
      if (type %in% c("unrooted", "radial") || !use.edge.length || 
          is.null(x$root.edge) || !x$root.edge) 
          root.edge <- FALSE
      phyloORclado <- type %in% c("phylogram", "cladogram")
      horizontal <- direction %in% c("rightwards", "leftwards")
      xe <- x$edge
      if (phyloORclado) {
          phyOrder <- attr(x, "order")
          if (is.null(phyOrder) || phyOrder != "cladewise") {
              x <- ape:::reorder.phylo(x)
              if (!identical(x$edge, xe)) {
                  ereorder <- match(x$edge[, 2], xe[, 2])
                  if (length(edge.color) > 1) {
                    edge.color <- rep(edge.color, length.out = Nedge)
                    edge.color <- edge.color[ereorder]
                  }
                  if (length(edge.width) > 1) {
                    edge.width <- rep(edge.width, length.out = Nedge)
                    edge.width <- edge.width[ereorder]
                  }
                  if (length(edge.lty) > 1) {
                    edge.lty <- rep(edge.lty, length.out = Nedge)
                    edge.lty <- edge.lty[ereorder]
                  }
              }
          }
          yy <- numeric(Ntip + Nnode)
          TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
          yy[TIPS] <- 1:Ntip
      }
      z <- ape:::reorder.phylo(x, order = "postorder") ##
      if (phyloORclado) {
          if (is.null(node.pos)) 
              node.pos <- if (type == "cladogram" && !use.edge.length) 
                  2
              else 1
          if (node.pos == 1) 
              yy <- .nodeHeight(z$edge, Nedge, yy)
          else {
              ans <- .C(ape::node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
                  1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
                  double(Ntip + Nnode), as.double(yy))
              xx <- ans[[5]] - 1
              yy <- ans[[6]]
          }
          if (!use.edge.length) {
              if (node.pos != 2) 
                  xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                    node.depth) - 1
              xx <- max(xx) - xx
          }
          else {
              xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                  z$edge.length)
          }
      }
      else {
          twopi <- 2 * pi
          rotate.tree <- twopi * rotate.tree/360
          if (type != "unrooted") {
              TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
              xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                  length.out = Ntip)
              theta <- double(Ntip)
              theta[TIPS] <- xx
              theta <- c(theta, numeric(Nnode))
          }
          switch(type, fan = {
              theta <- .nodeHeight(z$edge, Nedge, theta)
              if (use.edge.length) {
                  r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                    Nedge, z$edge.length)
              } else {
                  r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
                  r <- 1/r
              }
              theta <- theta + rotate.tree
              if (root.edge) r <- r + x$root.edge
              xx <- r * cos(theta)
              yy <- r * sin(theta)
          }, unrooted = {
              nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
              XY <- if (use.edge.length) ape::unrooted.xy(Ntip, Nnode, 
                  z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                  Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
              xx <- XY$M[, 1] - min(XY$M[, 1])
              yy <- XY$M[, 2] - min(XY$M[, 2])
          }, radial = {
              r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
              r[r == 1] <- 0
              r <- 1 - r/Ntip
              theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
              xx <- r * cos(theta)
              yy <- r * sin(theta)
          })
      }
      if (phyloORclado) {
          if (!horizontal) {
              tmp <- yy
              yy <- xx
              xx <- tmp - min(tmp) + 1
          }
          if (root.edge) {
              if (direction == "rightwards") 
                  xx <- xx + x$root.edge
              if (direction == "upwards") 
                  yy <- yy + x$root.edge
          }
      }
      if (no.margin) 
          par(mai = rep(0, 4))
      if (show.tip.label) 
          nchar.tip.label <- nchar(x$tip.label)
      max.yy <- max(yy)
      getLimit <- function(x, lab, sin, cex) {
          s <- strwidth(lab, "inches", cex = cex)
          if (any(s > sin)) 
              return(1.5 * max(x))
          Limit <- 0
          while (any(x > Limit)) {
              i <- which.max(x)
              alp <- x[i]/(sin - s[i])
              Limit <- x[i] + alp * s[i]
              x <- x + alp * s
          }
          Limit
      }
      if (is.null(x.lim)) {
          if (phyloORclado) {
              if (horizontal) {
                  xx.tips <- xx[1:Ntip]
                  if (show.tip.label) {
                    pin1 <- par("pin")[1]
                    tmp <- getLimit(xx.tips, x$tip.label, pin1, 
                      cex)
                    tmp <- tmp + label.offset
                  }
                  else tmp <- max(xx.tips)
                  x.lim <- c(0, tmp)
              }
              else x.lim <- c(1, Ntip)
          }
          else switch(type, fan = {
              if (show.tip.label) {
                  offset <- max(nchar.tip.label * 0.018 * max.yy * 
                    cex)
                  x.lim <- range(xx) + c(-offset, offset)
              } else x.lim <- range(xx)
          }, unrooted = {
              if (show.tip.label) {
                  offset <- max(nchar.tip.label * 0.018 * max.yy * 
                    cex)
                  x.lim <- c(0 - offset, max(xx) + offset)
              } else x.lim <- c(0, max(xx))
          }, radial = {
              if (show.tip.label) {
                  offset <- max(nchar.tip.label * 0.03 * cex)
                  x.lim <- c(-1 - offset, 1 + offset)
              } else x.lim <- c(-1, 1)
          })
      }
      else if (length(x.lim) == 1) {
          x.lim <- c(0, x.lim)
          if (phyloORclado && !horizontal) 
              x.lim[1] <- 1
          if (type %in% c("fan", "unrooted") && show.tip.label) 
              x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                  cex)
          if (type == "radial") 
              x.lim[1] <- if (show.tip.label) 
                  -1 - max(nchar.tip.label * 0.03 * cex)
              else -1
      }
      if (phyloORclado && direction == "leftwards") 
          xx <- x.lim[2] - xx
      if (is.null(y.lim)) {
          if (phyloORclado) {
              if (horizontal) 
                  y.lim <- c(1, Ntip)
              else {
                  pin2 <- par("pin")[2]
                  yy.tips <- yy[1:Ntip]
                  if (show.tip.label) {
                    tmp <- getLimit(yy.tips, x$tip.label, pin2, 
                      cex)
                    tmp <- tmp + label.offset
                  }
                  else tmp <- max(yy.tips)
                  y.lim <- c(0, tmp)
              }
          }
          else switch(type, fan = {
              if (show.tip.label) {
                  offset <- max(nchar.tip.label * 0.018 * max.yy * 
                    cex)
                  y.lim <- c(min(yy) - offset, max.yy + offset)
              } else y.lim <- c(min(yy), max.yy)
          }, unrooted = {
              if (show.tip.label) {
                  offset <- max(nchar.tip.label * 0.018 * max.yy * 
                    cex)
                  y.lim <- c(0 - offset, max.yy + offset)
              } else y.lim <- c(0, max.yy)
          }, radial = {
              if (show.tip.label) {
                  offset <- max(nchar.tip.label * 0.03 * cex)
                  y.lim <- c(-1 - offset, 1 + offset)
              } else y.lim <- c(-1, 1)
          })
      }
      else if (length(y.lim) == 1) {
          y.lim <- c(0, y.lim)
          if (phyloORclado && horizontal) 
              y.lim[1] <- 1
          if (type %in% c("fan", "unrooted") && show.tip.label) 
              y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                  cex)
          if (type == "radial") 
              y.lim[1] <- if (show.tip.label) 
                  -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
              else -1
      }
      if (phyloORclado && direction == "downwards") 
          yy <- y.lim[2] - yy
      if (phyloORclado && root.edge) {
          if (direction == "leftwards") 
              x.lim[2] <- x.lim[2] + x$root.edge
          if (direction == "downwards") 
              y.lim[2] <- y.lim[2] + x$root.edge
      }
      asp <- if (type %in% c("fan", "radial", "unrooted")) 
          1
      else NA
      L <- list(type = type, use.edge.length = use.edge.length, 
          node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
          show.node.label = show.node.label, font = font, cex = cex, 
          adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
          x.lim = x.lim, y.lim = y.lim, direction = direction, 
          tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
          align.tip.label = align.tip.label)
      invisible(L)
  }