#' Cluster subsamples of the data
#'
#' Given a nxp data matrix, this funciton will subsample the rows (samples),
#' cluster the subsamples, and return a nxn matrix with the probability of
#' co-occurance
#'
#' @param x the data on which to run the clustering (samples in rows)
#' @param k number of clusters to find for each clustering of a subsample (passed to clusterFunction)
#' @param clusterFunction a function that clusters a nxp matrix of data. Can also be given character values 'pam' or 'kmeans' to indicate use of internal wrapper functions. Must accept arguments 'x' and 'k' (whether uses them or not). See Details for format of what must return.
#' @param clusterArgs a list of parameter arguments to be passed to clusterFunction.
#' @param resamp.num the number of subsamples to draw
#' @param samp.p the proportion of samples to sample for each subsample.
#' @param classifyMethod method for determining which samples should be used in  the co-occurance matrix. "All"= all samples, "OutOfSample"= those not subsampled, and "InSample"=those in the subsample.  "All" and "OutOfSample" require that you provide classifyFunction to define how to classify those samples not in the subsample into a cluster. If "All" is chosen, all samples will be classified into clusters via the classifyFunctions, not just those that are out-of-sample. Note if not choose 'All' possible to get NAs in resulting D matrix (particularly if not enough subsamples taken).
#' @param classifyFunction a function which, given the output of clusterFunction and new data points, will classify the new data points into a cluster.
#' @param ncores integer giving the number of cores. If ncores>1, mclapply will be called.
#' @param ... arguments passed to mclapply (if ncores>1)
#'
#' @details The clusterFunction must be a function that takes as an argument 'x' which is a nxp matrix  of data and integer 'k'. It minimally must return a list with element named 'clustering' giving the vector of cluster ids. To be incorporated with the larger hierarchy, it should be list with elements of a partition object, just as is returned by pam. Generally, the user will need to write a wrapper function to do this. In the case of pam or kmeans, the user can identify clusterFunction as "pam" or "kmeans", and the package functions will use internally written wrappers for the clusterFunction and classifyFunction arguments. Additional arguments should be supplied via clusterArgs.
#'
#' @details The classifyFunction should take as an object a data matrix 'x' with samples on the rows, and the output of the clusterFunction. Note that the function should assume that the input 'x' is not the same samples that were input to the clusterFunction (but can assume that it is the same number of features/columns)
#'
#' @return A nxn matrix of co-occurances
#'
#' @examples
#' data(simData)
#' subD<-subsampleClustering(simData,k=3,clusterFunction="kmeans",
#' clusterArgs=list(nstart=10),resamp.n=100,samp.p=0.7)
#' heatmap(subD)
#' @export
subsampleClustering<-function(x,k,clusterFunction="pam", clusterArgs=NULL, classifyMethod=c("All","InSample","OutOfSample"),classifyFunction=NULL,resamp.num = 100, samp.p = 0.7,ncores=1,... )
{
	if(!is.function(clusterFunction)){
		if(clusterFunction%in%c("pam","kmeans")){
			if(clusterFunction=="pam"){
				clusterFunction<-.pamCluster
				classifyFunction<-.pamClassify
			}
			else if(clusterFunction=="kmeans"){
				clusterFunction<-.kmeansCluster
				classifyFunction<-.kmeansClassify
			}
		}
		else stop("clusterFunction must define a function")
	}
	classifyMethod<-match.arg(classifyMethod)
	if(classifyMethod %in% c("All","OutOfSample") && missing(classifyFunction)){
		stop("Must provide a classification function for the 'All' or 'OutOfSample' options")
	}
	N <- dim(x)[1]
    subSize <- round(samp.p * N)
    idx<-replicate(resamp.num,sample(1:N,size=subSize)) #each column a set of indices for the subsample.
	perSample<-function(ids){
		result<-do.call(clusterFunction,c(list(x=x[ids,],k=k),clusterArgs))
		if(classifyMethod=="All") classX<-classifyFunction(x,result)
		if(classifyMethod=="OutOfSample"){
			classElse<-classifyFunction(x[-ids,], result)
			classX<-rep(NA,N)
			classX[-ids]<-classElse
		}
		if(classifyMethod=="InSample"){
			classX<-rep(NA,N)
			classX[ids]<-result$clustering
		}
        D <- outer(classX, classX, function(a, b) a == b)
		Dinclude<-matrix(1,N,N)
		whNA<-which(is.na(classX))
		if(length(whNA)>0){
			Dinclude[whNA,]<-0 #don't add them to the denominator either
			Dinclude[,whNA]<-0
			D[whNA,]<-0 #don't add to sum
			D[,whNA]<-0
		}
		return(list(D=D,Dinclude=Dinclude))
	}
	if(ncores==1){
		DList<-apply(idx,2,perSample)
	}
	else{
		DList<-parallel::mclapply(1:ncol(idx),function(nc){perSample(idx[,nc])},mc.cores=ncores,...)
	}
	DDenom<-Reduce("+",lapply(DList,function(y){y$Dinclude}))
	DNum<-Reduce("+",lapply(DList,function(y){y$D}))
	Dbar = DNum/DDenom
	rownames(Dbar)<-colnames(Dbar)<-rownames(x)
	return(Dbar)
}
