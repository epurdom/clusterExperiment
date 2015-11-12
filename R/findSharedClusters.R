#' Find sets of samples that stay together across clusterings
#' Find sets of samples that stay together across clusterings in order to define a new clustering vector.
#'  
#' @param clusterMat A matrix with samples on rows and different clustering assignments on each column.
#' @param clusterFunction the clustering to use (passed to \code{\link{clusterD}}); currently must be of type '01'. 
#' @param minSize minimum size required for a set of samples to be considered in a cluster because of shared clustering, passed to \code{\link{clusterD}}
#' @param proportion The proportion of times that two sets of samples should be together in order to be grouped into a cluster
#' @plot logical. If true, will create a heatmap of the shared percentages with the final cluster found by function
#'
#' @return A list with values 
#' \itemize{
#'
#' \item{\code{clustering}}{vector of cluster assignments, with "-1" implying unassigned}
#'
#' \item{\code{percentageShared}}{a nsample x nsample matrix of the percent co-occurance across clusters used to find the final clusters}
#' }


findSharedClusters<-function(clusterMat,proportion=1,clusterFunction="hierarchical",plot=FALSE,minSize=5){
	#clusterMat a nsample x nclusteringMethod matrix of integered valued clusterings
	#minSize determines minimimum size needed to declare cluster replicated
	singleValueClusters<-apply(clusterMat,1,paste,collapse=";")
	allUnass<-paste(rep("-1",length=ncol(clusterMat)),collapse=";")
	if(proportion==1){
		uniqueSingleValueClusters<-unique(singleValueClusters)
		tab<-	table(singleValueClusters)
		tab<-tab[tab>=minSize]
		tab<-tab[names(tab)!=allUnass]
		cl<-match(singleValueClusters,names(tab))
		cl[is.na(cl)]<- -1		
		return(cl)
	}
	else{
		if(is.character(clusterFunction)){
			typeAlg<-.checkAlgType(clusterFunction)
			if(typeAlg!="01") stop("findSharedClusters is only implemented for '01' type clustering functions (see help of clusterD)")
		}
		sharedPerct<-diag(rep(1,length=nrow(clusterMat)))
		#is this any faster?
		if(require(gtools)){
			clusterMat<-data.frame(t(clusterMat))
			indMat<-gtools::combinations(n=ncol(clusterMat),r=2)
			tmp<-unlist(mapply(clusterMat[indMat[,1]],clusterMat[indMat[,2]], FUN=function(x,y){
				v<-as.character(x)!="-1" & as.character(y)!="-1"
				w<-x==y & v
				return(sum(w)/sum(v))
				},SIMPLIFY=FALSE))
			tmp[is.na(tmp)]<-0
			sharedPerct[lower.tri(sharedPerct)]<-tmp
			sharedPerct<-t(sharedPerct)
			sharedPerct[lower.tri(sharedPerct)]<-tmp
			sharedPerct[sharedPerct<proportion]<-0
			cl<-clusterD(D=sharedPerct,clusterFunction=clusterFunction,alpha = 1-proportion, minSize=minSize, format="vector")
			# clB<-findBlocks(Dbar=sharedPerct,alpha = 1-proportion, minSize=minSize, minSize.core=2, format="vector",coreValue=proportion)
			# clB2<-findBlocks(Dbar=sharedPerct,alpha = 1-proportion, findCoreType="track",minSize=minSize, minSize.core=2, format="vector",coreValue=proportion)
			# browser()
			# par(mfrow=c(1,2))
			if(plot && require(NMF)) NMF::aheatmap(sharedPerct,annCol=data.frame("Cluster"=factor(cl)),Colv="Rowv",annColors=list(bigPalette))
			# aheatmap(sharedPerct,annCol=data.frame("Cluster"=factor(clB2)),Colv="Rowv",annColors=list(brainUtils:::.thisPal))
			
		}
		
		else{ stop("Must have gtools installed")}
		# for(ii in 1:(ncol(sharedPerct)-1)){
		# 	for(jj in (ii+1):ncol(sharedPerct)){
		# 		x<-clusterMat[ii,]
		# 		y<-clusterMat[jj,]
		# 		share<-which(as.character(x)!="-1" & as.character(y)!="-1")
		# 		if(length(share)==0){
		# 			sharedPerct[ii,jj]<-sharedPerct[jj,ii]<-0
		# 		}
		# 		else{
		# 			x<-x[share]
		# 			y<-y[share]
		# 			sharedPerct[ii,jj]<-sharedPerct[jj,ii]<-length(which(x==y))/length(y)
		# 		}
		# 	}
		# }
		##Now create clusters via...??
		return(list(clustering=cl,percentageShared=sharedPerct))
	}
}