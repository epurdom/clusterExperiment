### input to clustering:
# pam : x or dis
# hier : dis
# kmeans : x
# spectral (SamSPECTRAL for flow cytometry; kernlab for standard; kknn for similarity based on knn rather than kmeans): kernlab is either x or a kernel function
###what rdname should this be? Not S4 methods...
#' @rdname ClusterFunction-methods
#' @export
.builtInClusterFunctions<-c("pam"=pamObj,"kmeans"=kmeansObj)
builtInClusterFunctions<-names(.builtInClusterFunctions)
################
##Internal wrapper functions for kmeans and pam
################
kmeansObj<-clusterFunction(clusterFUN=.kmeansCluster, classifyFUN=.kmeansClassify,inputType="X",inputClassifyType="X",algorithmType="K")

###Kmeans
.kmeansCluster <- function(x,k, checkArgs,cluster.only,...) { 
    passedArgs<-.getPassedArgs(FUN=stats::kmeans,passedArgs=list(...),checkArgs=checkArgs)
	  out<-do.call(stats::kmeans,c(list(x=t(x),centers=k),passedArgs)
  if(cluster.only) return(out$cluster)
  else return(.kmeansPartitionObject(x,out)) 
} 
.kmeansClassify <- function(x, clusterResult) { 
  centers <- clusterResult$mediods
  suppressWarnings(stats::kmeans(t(x), centers, iter.max = 1, algorithm = "Lloyd")$cluster) #probably uses this so always classifies points to centers
} 
#make partition object same form as pam output
.kmeansPartitionObject<-function(x,kmeansObj){ 
  dissE<-(cluster::daisy(t(x)))^2
  silObj<-cluster::silhouette(kmeansObj$cl,dissE^2)
  silinfo<-list(widths=silObj, clus.avg.widths=summary(silObj)$clus.avg.widths, ave.width=summary(silObj)$avg.width)
  return(list(mediods=kmeansObj$centers, clustering=kmeansObj$cluster, call=NA,silinfo=silinfo, objective=NA, diss=dissE, data=x))
}

###Pam
pamObj<-clusterFunction(clusterFUN=.pamCluster, classifyFUN=.pamClassify,inputType="either",inputClassifyType="X",algorithmType="K")

.pamCluster<-function(x,diss,k,checkArgs,cluster.only,...){
      passedArgs<-.getPassedArgs(FUN=cluster::pam,passedArgs=list(...),checkArgs=checkArgs)
	  input<-.checkXDissInput(x,diss,checkDiss=FALSE){
	  if(input=="X") return(do.call(cluster::pam, c(list(x=x,k=k, cluster.only=cluster.only), passedArgs)))
      if(input=="diss" | input=="both") return(do.call(cluster::pam, c(list(x=D,k=k, diss=TRUE, cluster.only=cluster.only), passedArgs)))
    }
.pamClassify <- function(x, clusterResult) { #x p x n matrix
  .genericClassify(x,clusterResult$medoids)
} 
.genericClassify<-function(x,centers){
    innerProd<-tcrossprod(t(x),centers) #a n x k matrix of inner-products between them
    distMat<-as.matrix(dist(rbind(t(x),centers)))
    distMat<-distMat[1:ncol(x),(ncol(x)+1):ncol(distMat)]
    apply(distMat,1,which.min)	
}
# .hierCluster<-function(x,k,...){
# 	argList<-list(...)
# 	hout<-do.call("hclust",c(list(dist(x)),argList))
# 	stats::cutree(tree, k = k)
#
# }
# .hierClassify<-function(x,clusterResult){
#
# }


