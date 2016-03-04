#' @title Make hierarchy of set of clusters
#'
#' @description Makes a dendrogram of a set of clusters based on hclust on the mediods of the cluster.
#'
#'
#' @param dat data to define the mediods from. Samples on rows, features on columns
#' @param cl A numeric vector with cluster assignments to compare to clRef.``-1'' indicates the sample was not assigned to a cluster.
#' @param full Logical as to whether to expand to include all samples as tips (useful for heatmap), or to have dendrogram just of the clusters (useful for getBestGenes)
#' @param unassigned what to do with with unassigned ("-1") cluster; only relevant if full=TRUE. Default is to remove those samples, so that only get heatmap of clustered observations. ``outgroup'' puts the unassigned into an outgroup on the dendrogram. ``cluster'' clusters them with the mediods. 
#' @param ... arguments passed to hclust
#' @details If full=TRUE, then the expanded dendrogram is created by giving all of the samples the mediod of the cluster and applying hclust to it; if further unassigned=="cluster", then the expanded dendrogram is created by hclust of the expanded mediod data plus the original unclustered observations.
#'
#'
#' @examples
#' data(simData)
#' #create a clustering, for 8 clusters (truth was 3)
#' cl<-clusterAll(simData,clusterFunction="pam",subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))$cl
#' #create dendrogram of clusters:
#' hcl<-clusterHclust(dat=simData,cl,full=FALSE)
#' plot(hcl)
#'
#' #create dendrogram for plotting with data in heatmap:
#' hclData<-clusterHclust(dat=simData,cl,full=TRUE)
#' dualHeatmap(cl,heatData=simCount,clusterData=hclData,colorScale=seqPal5,
#'	annCol=data.frame(PAM8=cl,Truth=trueCluster))
 
clusterHclust<-function(dat,cl,full=TRUE,unassigned=c("remove","outgroup","cluster"),...){
	unassigned<-match.arg(unassigned)
	if(nrow(dat)!=length(cl)) stop("cl must be the same length as the number of rows of dat")
	if(is.null(rownames(dat))) rownames(dat)<-as.character(1:nrow(dat))
	if(is.factor(cl)){warning("cl is a factor. Converting to numeric, which may not result in valid conversion")
			cl<-as.numeric(as.character(cl))
	}
	whRm<- which(cl!=-1)
	clFactor<-factor(cl[whRm])
	mediods<-do.call("rbind",by(dat[whRm,],clFactor,function(x){apply(x,2,median)}))
	rownames(mediods)<-levels(clFactor)
	nPerCluster<-table(clFactor)
	
	if(full){
		#make fake data with just mediods as value per sample:
		fakeData<-do.call("rbind",lapply(levels(clFactor),function(x){
			ind<-which(clFactor==x) #indices of those in this cluster
			med<-mediods[x,]
			mat<-matrix(rep(med,length(ind)),nrow=length(ind),byrow=TRUE)
			rownames(mat)<-rownames(dat[whRm,])[ind]
			return(mat)
		}))
		fakeData<-fakeData[rownames(dat[whRm,]),]
		if(length(whRm)!=nrow(dat) && unassigned != "remove"){
			if(unassigned=="outgroup"){
				#hard to use merge and then get the indices to go back to the same ones
				#cheat and add large amount to the unassigned so that they don't cluster to 
				outlierDat<-dat[-whRm,]
				maxAss<-max(dat[whRm,])
				outlierDat<-outlierDat+maxAss+10e6
				fakeData<-rbind(fakeData,outlierDat)
				fakeData<-fakeData[rownames(dat),]
				# fullD<-as.dendrogram(stats::hclust(dist(fakeData)))
				# unassD<-as.dendrogram(stats::hclust(dist(dat[-whRm,])))
				# return(merge(fullD,unassD))

			}
			if(unassigned=="cluster"){
				#add remaining to fake data and let them cluster
				fakeData<-rbind(fakeData,dat[-whRm,])
				fakeData<-fakeData[rownames(dat),]
				return(as.dendrogram(stats::hclust(dist(fakeData))))
			}
		}
		fullD<-as.dendrogram(stats::hclust(dist(fakeData)^2),...)
		if(length(whRm)!=nrow(dat) && unassigned == "outgroup"){
			#need to get rid of super long outgroup arm
			
			armLength<-max(attributes(fullD[[1]])$height,attributes(fullD[[2]])$height)
			attributes(fullD)$height<-armLength+.1*armLength
		}
		return(fullD)
	}
	else{
		return(as.dendrogram(stats::hclust(dist(mediods)^2,members=nPerCluster,...)))		
	}
}

#' @title Merge clusters based on dendrogram
#'
#' @description Takes an input of hierarchical clusterings of clusters and returns estimates of number of proportion of non-null and merges those below a certain cutoff.
#'
#'
#' @param dat data to perform the test on. Samples on rows, features on columns
#' @param cl A numeric vector with cluster assignments to compare to clRef.``-1'' indicates the sample was not assigned to a cluster.
#' @param dendro If provided, is dendrogram providing hierarchical clustering of clusters in cl; mainly useful to speed up calculations if clusterHclust has already been called on the (cl,dat) pair. The default (NULL) is to recalculate it with the given (cl,dat) pair.
#' @param mergeMethod method for calculating proportion of non-null that will be used to merge clusters (if 'none', no merging will be done). See details for description of methods. 
#' @param cutoff cutoff for merging clusters based on the proportion of non-nulls in the comparison of the clusters (i.e. value between 0,1, where lower values will make it harder to merge clusters).
#' @param plotType what type of plotting of dendrogram. If 'all', then all the estimates of proportion non-null will be plotted; if 'mergeMethod', then only the value used in the merging is plotted for each node. 
#' @param ... arguments passed to the \code{\link{plot.phylo}} function of \code{ade4} that plots the dendrogram. 
#' 
#' @details "JC" refers to the method of Ji and Cai (2007), and implementation of "JC" method is copied from code available on Jiashin Ji's website, December 16, 2015 (http://www.stat.cmu.edu/~jiashun/Research/software/NullandProp/). "locfdr" refers to the method of Efron (2004) and is implemented in the package \code{\link{locfdr}}. "MB" refers to the method of Meinshausen and Buhlmann (2005) and is implemented in the package \code{\link{howmany}}. "adjP" refers to the proportion of genes that are found significant based on a FDR adjusted p-values (method "BH") and a cutoff of 0.05. 
#'
#' @details If \code{mergeMethod} is not equal to 'none' then the plotting will indicate the merged clusters (assuming \code{plotType} is not 'none'). 
#' @examples
#' data(simData)
#'  #create a clustering, for 8 clusters (truth was 3)
#'  cl<-clusterAll(simData,clusterFunction="pam",subsample=FALSE,
#'  sequential=FALSE, clusterDArgs=list(k=8))$cl
#' #merge clusters with plotting. Note argument 'use.edge.length' can improve readability
#' mergeResults<-mergeClusters(simData,cl=cl,plot=TRUE,plotType="all",
#' mergeMethod="adjP",use.edge.length=FALSE)
#' #compare merged to original on heatmap
#' hclData<-clusterHclust(dat=simData,cl,full=TRUE)
#' dualHeatmap(cl,heatData=simCount,clusterData=hclData,colorScale=seqPal5,
#'	annCol=data.frame(Original=cl,Merged=mergeResults$cl,Truth=trueCluster))

 
mergeClusters<-function(dat,cl,dendro=NULL,mergeMethod=c("none","adjP","locfdr","MB","JC"),cutoff=0.1,plotType=c("none","all","mergeMethod"),countData=TRUE,...){
	if(is.factor(cl)){warning("cl is a factor. Converting to numeric, which may not result in valid conversion")
			cl<-as.numeric(as.character(cl))
	}
	if(!is.null(dendro)){
		#check valid
		ncluster<-length(table(cl[cl>0]))
		if(nobs(dendro)!=ncluster){
			warning("Not a valid input dendrogram (not equal to the number of non -1 clusters in cl). Will recompute dendrogram")
			dendro<-NULL
		}
	}
	if(is.null(dendro)){
		dendro<-if(countData) clusterHclust(dat=log(dat+1),cl,full=FALSE) else clusterHclust(dat=dat,cl,full=FALSE)		
	}
  	mergeMethod<-match.arg(mergeMethod)
  	plotType<-match.arg(plotType)
	if(plotType=="mergeValue" & mergeMethod=="none") stop("can only plot merge method values if one method is selected")
	#get test-statistics for the contrasts corresponding to each node (and return all)
	sigTable<-getBestGenes(cl,dat,type=c("Dendro"),dendro=dendro,returnType=c("Table"),contrastAdj=c("All"),number=ncol(dat),p.value=1,voomCorrection=countData)
	#divide table into each node.
	sigByNode<-by(sigTable,sigTable$ContrastName,function(x){
		mb<-.myTryFunc(pvalues=x$P.Value,FUN=.m1_MB)
		locfdr<-.myTryFunc(tstats=x$t,FUN=.m1_locfdr)
		jc<-.myTryFunc(tstats=x$t,FUN=.m1_JC)
		return(c("adjP"=.m1_adjP(x$adj),"locfdr"=locfdr,"MB"=mb,"JC"=jc))
	})
	newcl<-cl
	phylo4Obj<-.makePhylobaseTree(dendro,"dendro")
	
	if(mergeMethod!="none"){
		#go up tree and merge clusters
		valsPerNode<-sapply(sigByNode,function(x){signif(x[[mergeMethod]],2)})
		nodesBelowCutoff<-names(valsPerNode)[which(valsPerNode<cutoff)] #names of nodes below cutoff
		
		#find nodes where *all* descendants are below cutoff
		allTipNames<-phylobase::labels(phylo4Obj)[phylobase::getNode(phylo4Obj,  type=c("tip"))]
		whToMerge<-sapply(nodesBelowCutoff,function(node){
			desc<-phylobase::descendants(phylo4Obj,node,type = c("all"))
			return(all(names(desc) %in% nodesBelowCutoff | names(desc) %in% allTipNames))
		})
		if(length(whToMerge)>0){
			nodesToMerge<-nodesBelowCutoff[whToMerge]

			#now find top ones
			whAnc<-sapply(nodesToMerge,function(node){
				anc<-phylobase::ancestors(phylo4Obj,node,type="all")
				return(!any(names(anc) %in% nodesToMerge))
			})
			nodesAtTop<-nodesToMerge[whAnc]
		
			#make new clusters
			temp<-lapply(nodesAtTop,function(node){
				tips<-phylobase::descendants(phylo4Obj,node,type="tips")
				if(any(!names(tips) %in% as.character(cl))) stop("coding error-- tips don't match values of cl")
					newcl[cl%in% tips]<<- as.numeric(tips[[1]])
			})
			#make consecutive integers
			newcl<-as.numeric(factor(newcl,levels=unique(cl[cl>0])))
			#deal with -1/-2
			newcl[is.na(newcl)]<-cl[is.na(newcl)]
		}
	}
	if(plotType!="none"){
		# phylobase has bug in plotting! submitted to their github
		# move to ape package...	
		phyloObj<-as(phylo4Obj, "phylo")

		#####
		#convert names of internal nodes for plotting
		#####
		allInternal<-phyloObj$node
		#match to order of tree
		m<-match(allInternal,names(sigByNode))
		edgeLty<-rep(1,nrow(phyloObj$edge))
		if(mergeMethod!="none" & length(whToMerge)>0){
			#browser()
			whMerge<-which(phyloObj$node.label %in% nodesToMerge)#which of nodes merged
			nodeNumbers<-(length(phyloObj$tip)+1):max(phyloObj$edge)
			whEdge<-which(phyloObj$edge[,1] %in% nodeNumbers[whMerge])
			edgeLty[whEdge]<-2
		}
		if(plotType=="mergeMethod"){
				phyloObj$node.label<-as.character(valsPerNode)
		}
		if(plotType=="all") phyloObj$node.label<-sapply(sigByNode[m],function(x){paste(paste(names(x),signif(x,2),sep=":"),collapse=",\n")})
		ape::plot.phylo(phyloObj,show.node=TRUE,edge.lty=edgeLty,...)
	}
	nodePropTable<-do.call("rbind",sigByNode)
	nodePropTable<-data.frame("Node"=names(sigByNode),"Contrast"=sigTable$Contrast[match(names(sigByNode),sigTable$ContrastName)],nodePropTable)

	invisible(list(cl=newcl,oldClToNew=table(Original=cl,New=newcl),propDE=nodePropTable,originalClusterDendro=dendro))
}

.myTryFunc<-function(FUN,...){
	x<-try(FUN(...))
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
	musigma<-try(.EstNull.func(tstats))
	.epsest.func(tstats,musigma$mu,musigma$s) #gives proportion of non-null
}
