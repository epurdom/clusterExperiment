################
##Internal wrapper functions for kmeans and pam
################

###Kmeans
.kmeansCluster <- function(x,k, ...) { 
  out<-stats::kmeans(t(x),centers=k,...)
  out<-.kmeansPartitionObject(x,out) #make it a partition object like pam.
  #out$clustering<-out$cluster #stupid difference in naming...
  return(out)
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
  return(list(mediods=kmeansObj$centers,clustering=kmeansObj$cluster,call=NA,silinfo=silinfo,objective=NA,diss=dissE,data=x))
}

###Pam
.pamCluster <- function(x,k, ...) { cluster::pam(x=t(x),k=k,...) }  #x p x n matrix
.pamClassify <- function(x, clusterResult) { #x p x n matrix
  center<-clusterResult$medoids
  innerProd<-tcrossprod(t(x),center) #a n x k matrix of inner-products between them
  distMat<-as.matrix(dist(rbind(t(x),center)))
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


