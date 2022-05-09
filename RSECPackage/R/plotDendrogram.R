#'@title Plot dendrogram of RSECClass object
#'  
#'@description Extends plotDendrogram to RSECClass object
#'  
#'@param x a \code{\link{RSECClass}} object.
#'@param mergeInfo What kind of information about merge to plot on dendrogram. 
#'  If not equal to "none", will replicate the kind of plot that 
#'  \code{\link{mergeClusters}} creates, and the input to \code{mergeInfo} 
#'  corresponds to that of \code{plotInfo} in \code{mergeClusters}.
#' @return A dendrogram is plotted. Returns (invisibly) a list with elements
#' \itemize{
#' \item{\code{plottedObject}}{ the \code{phylo} object that is plotted.}
#' \item{\code{originalObject}}{ the \code{phylo} object before adjusting the
#' node/tip labels. }
#' }
#' @aliases plotDendrogram
#' @importFrom ape plot.phylo nodelabels
#' @seealso
#'   \code{\link{mergeClusters}},\code{\link[clusterExperiment]{plotDendrogram}}
#'   
#' @export
#' 
#' @examples
#' data(simData)
#' 
#' #create a clustering, for 8 clusters (truth was 3) 
#' cl <-clusterSingle(simData, subsample=FALSE, 
#' sequential=FALSE, 
#' mainClusterArgs=list(clusterFunction="pam", clusterArgs=list(k=8)))
#' 
#' #create dendrogram of clusters and then 
#' # merge clusters based ondendrogram: 
#' cl <- makeDendrogram(cl) 
#' cl <- mergeClusters(cl,mergeMethod="adjP",DEMethod="limma",
#'    cutoff=0.1,plot=FALSE) 
#' plotDendrogram(cl) 
#' plotDendrogram(cl,leafType="samples",whichClusters="all",plotType="colorblock")
#' 
#' @export
#' @rdname plotDendrogram
setMethod(
  f = "plotDendrogram",
  signature = "RSECClass",
  definition = function(x,mergeInfo="none", ...)
  {    
    possibleMergeValues<-c("none", "all","mergeMethod",.availMergeMethods)
    if(!is.null(x@merge_nodeProp)){
      otherVals<-colnames(x@merge_nodeProp)[!colnames(x@merge_nodeProp)
          %in% c("NodeId","Contrast")]
      possibleMergeValues<-unique(c(possibleMergeValues,otherVals))
      
    }
    mergeInfo<-match.arg(mergeInfo,possibleMergeValues)    
    if(is.na(x@merge_dendrocluster_index)) mergeInfo<-"none"
    if(mergeInfo=="none"){
      mergeInfo<-NULL
      mergeMethod<-NULL
    }
    else if(!is.na(x@merge_dendrocluster_index)){ 
      #has node prop info, even if no merge cluster
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
    mergeList<-list(
      mergeMethod=mergeMethod,
      mergePlotType=mergeInfo,
      mergeOutput=nodeMergeInfo(x))
    
    callNextMethod(x=x,mergeInfo=mergeList,...)
    
  })




########
# Internal plotting function used by both mergeClusters and plotDendrogram
#' @importFrom graphics plot
#' @importFrom ape plot.phylo phydataplot
.plotDendro<-function(dendro, leafType="clusters",mergePlotType=NULL,mergeMethod=NULL,mergeOutput=NULL,clusterLegendMat=NULL,clObj=NULL,plotType=c("name","colorblock"),removeOutbranch=FALSE,legend="below",clusterLabelAngle=45,...){
  plotType<-match.arg(plotType)
  outbranch<- "outbranch root" %in% phylobase::tdata(dendro)$Position
  
  ############
  ## Everything in this code assumes that phylobase::getNode(dend@dendro_samples,type="tip") is in same order as phyloObj$tip.label
  ############
    
  #---
  #remove the outbranch from the dendrogram and update clObj and mTipsToSamples
  #(note this is using phylo4 obj)
  #---	
  if(outbranch & removeOutbranch & leafType=="samples"){
	  ##Find node that is cluster child of root	  
    rootNode<-phylobase::rootNode(dendro)
    rootChild<-phylobase::descendants(dendro,node=rootNode,type="children")
	position<-.matchToDendroData(inputValue=rootChild,dendro=dendro,matchColumn="NodeIndex",returnColumn="Position")
	if(!any(position=="cluster hierarchy node")) stop("coding error -- child of root with outbranch isn't cluster hierarchy node")
	if(all(position=="cluster hierarchy node")) stop("coding error -- both child of root are 'cluster hierarchy node', but also have root is 'outbranchroot'")	
	clusterNode<-rootChild[which(position=="cluster hierarchy node")]
	
	#remove from clObj and update mTipsToSamples
	nSamples<-nTips(dendro) #save so have after changed dendro
    clusterTips<-phylobase::descendants(dendro,node=clusterNode,type="tip")
    if(length(clusterTips)==0) stop("Internal coding error: no unassigned samples in tree")
	whKeep<-.matchToDendroData(inputValue=clusterTips, dendro, matchColumn="NodeIndex", returnColumn="SampleIndex")
#	all(tipLabels(dendro)==row.names(clObj)[tdata(dendro,type="tip")$SampleIndex])

	####Subset dendro
	dendro<-phylobase::subset(dendro, node.subtree=clusterNode)
	ch<-.checkDendroSamplesFormat(dendro,checkLabels=FALSE)
	if(!is.logical(ch)) stop(ch)

	#need to create mTipsToSamples -- match of tips to samples ("SampleIndex" is index to full data)
	mTipsToSamplesOld <- .matchToDendroData(inputValue=phylobase::getNode(dendro,type="tip"), dendro, matchColumn="NodeIndex", returnColumn="SampleIndex")
	#these are still indices in the full sample clObj. Now need to get their indices in the subsetted one:
	mToSubset<-match(seq_len(nSamples),whKeep) #gives where each of old sample indices (1:n) map to in new order of future clObj (with NA for those that not in new clObj)
	mTipsToSamples<-mToSubset[mTipsToSamplesOld] #use that to map mTipsToSamplesNew to new order
	if(any(is.na(mTipsToSamples))) stop("coding error -- didn't update mTipsToSamples correctly")
		#these check against names, but doesn't always have names, etc.
	# if(!all(tipLabels(dendro)==row.names(clObj)[mTipsToSamplesOld])) stop("coding error -- names no longer match")
	# if(!all(sort(tipLabels(dendro))==sort(names(clObj[whKeep,])))) stop("coding error -- don't have the same set of tips after pruning outbranch")
	# if(!all(tipLabels(dendro)==row.names(clObj[whKeep,,drop=FALSE])[mTipsToSamples])) stop("coding error -- names no longer match")
    if(is.matrix(clObj)) clObj<-clObj[whKeep,,drop=FALSE] else clObj<-clObj[whKeep]		
    #set outbranch=FALSE because now doesn't exist in tree...
    outbranch<-FALSE
  }
  else{
	  if(leafType=="samples") mTipsToSamples <- .matchToDendroData(inputValue=phylobase::getNode(dendro,type="tip"), dendro, matchColumn="NodeIndex", returnColumn="SampleIndex")
	
  }
  
  #-------
  #convert to phylo object...
  #-------
  phyloObj <- .convertToPhyClasses(dendro, "phylo",convertNodes=TRUE,convertTips=(leafType!="samples")) 
  if(leafType=="samples"){
	  #not clear need this now that convert dendro it before send to .plotDendro
	  #these are in the order of dendro tips. To 
	  clNames<-phylobase::tipLabels(dendro)
  }
  plotArgs<-list(...)
  if(plotType=="colorblock" && is.null(clObj) && leafType=="samples") stop("Internal coding error: must provide a clustering if plotType='colorblock'")
  origPhylo<-phyloObj #so can return and get any information that gets overwritten in phyloObj
  
  #---------------
  ### For plotting of dendrogram for the merging
  ### Add information about the merging as node labels and change edge type
  ### Note: could probably have used nodelabels function and avoided some of this
  #---------------
  if(!is.null(mergeOutput)){
    annotNames<-c("NodeId","Contrast","isMerged", "mergeClusterId")
    methods<-colnames(mergeOutput)[!colnames(mergeOutput)%in%annotNames] #possible for which have proportion saved
  }
  doMerge<-!is.null(mergePlotType) && !is.null(mergeOutput) && mergePlotType %in% c("all",methods,"mergeMethod")
  if(doMerge){
    #####
    #convert names of internal nodes for plotting
    #####
    #match to order of tree
    whToMerge<-which(mergeOutput$isMerged)
    nodesToMerge<-as.character(mergeOutput$Node[whToMerge])
    m <- match( as.character(mergeOutput$Node),phyloObj$node.label)
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
  else{
	  phyloObj$node.label<-phylobase::nodeLabels(dendro) #use the user-given/assigned node labels
	
  }
  ###############
  ### Deal with clusterLegend object: 
  ### - Make default if not provided and 
  ### - If # of clusterings>1 make clusterLegend and clObj matrix appropriate
  ###############
  
  if(plotType=="colorblock"){
	  
    clusterLegend<-TRUE #doesn't do anything right now because phydataplot doesn't have option of no legend...
    if(is.null(clusterLegendMat)){ 
      #----
      #make default colors, works for vector or matrix clObj
      #----
      clusterIds<-sort(unique(as.vector(clObj)))
      clusterLegendMat <- cbind("clusterIds"=clusterIds, "name"=clusterIds, "color"=bigPalette[seq_along(clusterIds)])
    }
    else{
      if(is.matrix(clObj) && ncol(clObj)>1){
        #if not provide list of cluster legends, do only 1st clustering provided (temporary while fixing so works for matrix)
        if(!is.list(clusterLegendMat) ) clObj<-clObj[,1,drop=FALSE]
        else{
          #----
          #create one big clObj/clusterLegendMat object that will allow for coloring that is okay.
          #----
          nclusters<-ncol(clObj)
          if(length(clusterLegendMat)!=nclusters) stop("Internal coding error -- wrong length of colors for clustering")
          newClusterLegendMat<-clusterLegendMat[[1]]
          newCl<-clObj[,1]
          
          #make it general in case some day want more than just 2 clusterings
          for(ii in 2:nclusters){
            currMat<-clusterLegendMat[[ii]]
            currCl<-clObj[,ii]
            
            #note that because subset to those samples that are not -1/-2 on
            # clObj[,1], may have entire clusters in other columns of clObj that disappear but still in color matrix with no entry in clObj
            #reduce down the currMat to accomodate that
            whExist<-which(as.numeric(currMat[,"clusterIds"]) %in% currCl)
            currMat<-currMat[whExist, ,drop=FALSE]
            
            whExistingColor<-which(currMat[,"color"] %in% newClusterLegendMat[,"color"])
            updatedCurrCl<-currCl
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
              updatedCurrCl[!is.na(mexist)]<-newFullId[!is.na(mexist)]
              
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
              updatedCurrCl[!is.na(mexist2)]<-newFullId2[!is.na(mexist2)]
              
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
            #change them all here (before changed)
            currCl<-updatedCurrCl
            newCl<-cbind(newCl,currCl)
            
          }
          clusterLegendMat<-newClusterLegendMat
          colnames(newCl)<-colnames(clObj)
          rownames(newCl)<-rownames(clObj)
          clObj<-newCl
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
	  #don't need this now that convert before send to .plotDendro
      m<-match(gsub("ClusterId","",phyloObj$tip.label),clusterLegendMat[,"clusterIds"])
      if(any(is.na(m))) stop("clusterIds in clusterLegend do not match dendrogram labels")
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
  	  if(is.matrix(clObj) && ncol(clObj)>1){
        if(plotType=="colorblock"){
          colorMat<-apply(clObj,2,function(x){
            m<-match(x,clusterLegendMat[,"clusterIds"])
            clusterLegendMat[m,"name"]
          })
          if(any(dim(colorMat)!=dim(clObj))) stop("Internal coding error: dimensions of colorMat don't match input")
          dimnames(colorMat)<-dimnames(clObj)
          cols<-clusterLegendMat[,"color"]
          names(cols)<-clusterLegendMat[,"name"]
		  #match samples to order of the tree:
          colorMat<-colorMat[mTipsToSamples,]
        }
        tip.color<-"black"
      }
      else{
        if(is.matrix(clObj)) clObj<-clObj[,1]
        m<-match(clObj,clusterLegendMat[,"clusterIds"])
        tip.color<-clusterLegendMat[m,"color"]		
        if(plotType=="colorblock"){
		  colorMat<-matrix(clusterLegendMat[m,"name"],ncol=1)
		  colorMat<-colorMat[mTipsToSamples,,drop=FALSE]
          rownames(colorMat)<-clNames
          cols<-clusterLegendMat[,"color"]
          names(cols)<-clusterLegendMat[,"name"]
        }	
		else tip.color<-tip.color[mTipsToSamples]		
      }
      if(plotType=="colorblock"){
        #make only edges going to/from nodes in cluster hierarchy have edge.width>0
		#Note that some of origPhylo$node.label have NA for value; won't get unique value for these; but these aren't part of cluster hierarchy anyway (in fact could probably just pick those that are NA and use them, but this is safer(?))
		ntips<-length(phyloObj$tip.label)
		positionValue<-.matchToDendroData(inputValue=origPhylo$node.label, dendro=dendro, matchColumn="NodeId", returnColumn="Position")
		whClusterNode<-which(as.character(positionValue)%in% c("cluster hierarchy node","cluster hierarchy tip")) + ntips 
        whEdgePlot<-which(apply(phyloObj$edge,1,function(x){any(x %in% whClusterNode)}))
        edge.width<-rep(0,nrow(phyloObj$edge))
        edge.width[whEdgePlot]<-1
      }
    }
  }
  else tip.color<-"black"
	 
  #Couldn't do this before, because needed cluster id to safely match to clusterLegendMat. 
  if(leafType=="clusters") phyloObj$tip.label<-phylobase::tipLabels(dendro)

  #this next code is hack to deal with error sometimes get if very long edge length -- usually due to unusual distance, etc.
  # Divides edge lengths so not too large.
  if(max(phyloObj$edge.length)>1e6) phyloObj$edge.length <- phyloObj$edge.length / max(phyloObj$edge.length) 
  
  #---------------
  # PLOTTING
  #-------------
  #the percentage of total width will be dataPct/(1+dataPct+offsetPct)
  #tree=64%
  #data=32%
  #offset=3%
  dataPct<-0.5 
  offsetPct<-0.01 
  
  .pctCalculation<-function(data,offset){
	  denom<-1+data+offset
	  return(c("tree"=1/denom,data=data/denom,offset=offset/denom))
  }
  # > .pctCalculation(.5,.01)
#          tree        data      offset
#   0.662251656 0.331125828 0.006622517
  
  prohibitOptions<-c("tip.color","node.pos","edge.width","horizontal","direction","type")
  if(any(prohibitOptions %in% names(plotArgs))) stop("User cannot set following options to plot.phylo:",paste(prohibitOptions, collapse=","))
  plotArgs<-c(plotArgs,list(tip.color=tip.color,node.pos=2,edge.width=edge.width))	
  if(plotType=="name") do.call(ape::plot.phylo,c(list(phyloObj),plotArgs))
  else{
    #if colorblock
    if(length(grep("show.tip",names(plotArgs)))==0) plotArgs$show.tip.label<-FALSE
	sn<-grep("show.node",names(plotArgs))
	if(length(sn)>0 && plotArgs[[sn]] && !doMerge){
		offsetPct<-1/8	
		dataPct<-0.6
		# > .pctCalculation(.6,1/6)
# 		      tree       data     offset
# 		0.56603774 0.33962264 0.09433962
	}
			
	#just calculate needed width for tree...
    phyloPlotOut<-do.call(.calculatePlotPhylo,c(list(phyloObj,plot=FALSE),plotArgs))
    treeWidth<-phyloPlotOut$x.lim[2]
    #plot dendrogram:
    do.call(ape::plot.phylo,c(list(phyloObj,x.lim=treeWidth*(1+dataPct+offsetPct)),plotArgs))
    
    nclusters<-ncol(colorMat)
    colnames(colorMat)<-NULL		    
    colInput<-function(n){cols}
    width<-treeWidth*dataPct/(nclusters+.5)
	ape::phydataplot(x=colorMat, phy=phyloObj, style="mosaic",offset=treeWidth*offsetPct, width = width, border = NA, lwd = 3,legend = legend, funcol = colInput)
    
    if(nclusters>1 & !is.null(colnames(clObj))){
      xloc<-treeWidth+treeWidth*offsetPct+seq(from=0,by=width,length=nclusters)
      xloc<-xloc+width/2
      ypos<-par("usr")[4]-0.025*diff(par("usr")[3:4])	
      adj<-c(0,0)		
      if("cex" %in% names(list(...))) labcex<-list(...)[["cex"]]
      else labcex<-1	
      text(x=xloc,y=ypos,labels=colnames(clObj),srt=clusterLabelAngle,xpd=NA,adj=adj,cex=labcex)
      
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
      x <- ape::reorder.phylo(x)
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
    yy[TIPS] <- seq_len(Ntip)
  }
  z <- ape::reorder.phylo(x, order = "postorder") ##
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
        xx.tips <- xx[seq_len(Ntip)]
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
        yy.tips <- yy[seq_len(Ntip)]
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
  asp <- if (type %in% c("fan", "radial", "unrooted")) 1 else NA
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
            align.tip.label = align.tip.label)
  invisible(L)
}
