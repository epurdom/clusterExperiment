#' Find sets of samples that stay together across clusterings
#' Find sets of samples that stay together across clusterings in order to define a new clustering vector.
#'  
#' @param clusterMat A matrix with samples on rows and different clustering assignments on each column.
#' @param clusterFunction the clustering to use (passed to \code{\link{clusterD}}); currently must be of type '01'. 
#' @param minSize minimum size required for a set of samples to be considered in a cluster because of shared clustering, passed to \code{\link{clusterD}}
#' @param proportion The proportion of times that two sets of samples should be together in order to be grouped into a cluster (if <1, passed to clusterD via alpha=1-proportion)
#' @param propUnassigned samples with greater than this proportion of assignments equal to '-1' are assigned a '-1' cluster value.
#' @param plot logical. If true, will create a heatmap of the shared percentages with the final cluster found by function
#'
#' @details The function tries to find a consensus cluster across many different clusterings of the same sample. It does so by creating a nsample x nsample matrix of the percentage of co-occurance of each sample and then calling clusterD to cluster the co-occurance matrix. The function assumes that '-1' labels indicate clusters that are not assigned to a cluster. Co-occurance with the unassigned cluster is treated differently than other clusters. The percent co-occurance is taken only with respect to those clusterings where both samples were assigned. Then samples with more than \code{propUnassigned} values that are '-1' across all of the clusterings are assigned a '-1' regardless of their cluster assignment. 
#'
#' If the majority
#' @return A list with values 
#' \itemize{
#'
#' \item{\code{clustering}}{vector of cluster assignments, with "-1" implying unassigned}
#'
#' \item{\code{percentageShared}}{a nsample x nsample matrix of the percent co-occurance across clusters used to find the final clusters. Percentage is out of those not '-1'}
#' \item{\code{noUnassignedCorrection}{a vector of cluster assignments before samples were converted to '-1' because had >\code{propUnassigned} '-1' values (i.e. the direct output of the \code{clusterD} output.)}}
#' }


findSharedClusters<-function(clusterMat,proportion=1,clusterFunction="hierarchical",propUnassigned=.5,plot=FALSE,minSize=5){
	#clusterMat a nsample x nclusteringMethod matrix of integered valued clusterings
	#minSize determines minimimum size needed to declare cluster replicated
	if(proportion==1){
		singleValueClusters<-apply(clusterMat,1,paste,collapse=";")
		allUnass<-paste(rep("-1",length=ncol(clusterMat)),collapse=";")
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
		##Make clusterMat character, just in case
		clusterMat<-apply(clusterMat,2,as.character)
		clusterMat[clusterMat== "-1"]<-NA
		sharedPerct<-.hammingdist(t(clusterMat)) #works on columns. gives a nsample x nsample matrix back.
		sharedPerct[is.na(sharedPerct)|is.nan(sharedPerct)]<-0 #have no clusterings for which they are both not '-1'
		cl<-clusterD(D=sharedPerct,clusterFunction=clusterFunction,alpha = 1-proportion, minSize=minSize, format="vector",clusterArgs=list(evalClusterMethod=c("average")))
		if(plot && require(NMF)) NMF::aheatmap(sharedPerct,annCol=data.frame("Cluster"=factor(cl)),Colv="Rowv",annColors=list(bigPalette))
		
		if(is.character(cl)) stop("coding error -- clusterD should return numeric vector")
		##Now define as unassigned any samples with >= propUnassigned '-1' values in clusterMat
		whUnassigned<-which(apply(clusterMat,2,function(x){sum(x== -1)/length(x)>propUnassigned}))
		clUnassigned<-cl
		clUnassigned[whUnassigned]<- -1
		return(list(clustering=clUnassigned,percentageShared=sharedPerct,noUnassignedCorrection=cl))
	}
}


#from: https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/
.hammingdist <- function(X, Y,uniqValue=1539263) {
	#samples are in the rows!
    if ( missing(Y) ) {
        uniqs <- na.omit(unique(as.vector(X)))
		if(uniqValue %in% uniqs) stop("uniqValue (",uniqValue,") is in X")
		isobsX<-abs(is.na(X)-1)
		if(any(is.na(X))){
			X[is.na(X)]<-uniqValue
		}
    	U <- X == uniqs[1]
        H <- t(U) %*% U
        N <- t(isobsX) %*% (isobsX)
 		for ( uniq in uniqs[-1] ) {
            U <- X == uniq
            H <- H + t(U) %*% U
        }
    } else {
        uniqs <- na.omit(union(X, Y))
		if(uniqValue %in% uniqs) stop("uniqValue (",uniqValue,") is in either X or Y")
		isobsX<-abs(is.na(X)-1)
		if(any(is.na(X))){
			X[is.na(X)]<-uniqValue
		}
		isobsY<-abs(is.na(Y)-1)
		if(any(is.na(Y))){
			Y[is.na(Y)]<-uniqValue
		}
        H <- t(X == uniqs[1]) %*% (Y == uniqs[1])
        N <- t(isobsX) %*% (isobsY)
 		for ( uniq in uniqs[-1] ) {
			A<-t(X == uniq) %*% (Y == uniq)
           H <- H + A
        }
    }
    H/N
}

