#' @title Make hierarchy of set of clusters
#'
#' @description Makes a dendrogram of a set of clusters based on hclust on the mediods of the cluster.
#'
#'
#' @param dat data to define the mediods from. Samples on rows, features on columns
#' @param cl A vector with cluster assignments to compare to clRef.``-1'' indicates the sample was not assigned to a cluster.
#' @param full Logical as to whether to expand to include all samples as tips (useful for heatmap), or to have dendrogram just of the clusters (useful for getBestGenes)
#' @param unassigned what to do with with unassigned ("-1") cluster; only relevant if full=TRUE. Default is to remove those samples, so that only get heatmap of clustered observations. ``outgroup'' puts the unassigned into an outgroup on the dendrogram. ``cluster'' clusters them with the mediods. 
#'@details If full=TRUE, then the expanded dendrogram is created by giving all of the samples the mediod of the cluster and applying hclust to it; if further unassigned=="cluster", then the expanded dendrogram is created by hclust of the expanded mediod data plus the original unclustered observations.
#'
#'

clusterHclust<-function(dat,cl,full=TRUE,unassigned=c("remove","outgroup","cluster"),...){
	unassigned<-match.arg(unassigned)
	if(nrow(dat)!=length(cl)) stop("cl must be the same length as the number of rows of dat")
	if(is.null(rownames(dat))) rownames(dat)<-as.character(1:nrow(dat))
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
