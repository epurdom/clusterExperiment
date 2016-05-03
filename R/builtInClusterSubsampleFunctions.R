################
##Internal wrapper functions for kmeans and pam
################

###Kmeans
.kmeansCluster <- function(x,k, ...) { 
  out<-stats::kmeans(x,centers=k,...)
  out<-.kmeansPartitionObject(x,out) #make it a partition object like pam.
  #out$clustering<-out$cluster #stupid difference in naming...
  return(out)
} 
.kmeansClassify <- function(x, clusterResult) { 
  centers <- clusterResult$mediods
  suppressWarnings(stats::kmeans(x, centers, iter.max = 1, algorithm = "Lloyd")$cluster) #probably uses this so always classifies points to centers
} 
#make partition object same form as pam output
.kmeansPartitionObject<-function(x,kmeansObj){
  dissE<-(cluster::daisy(x))^2
  silObj<-cluster::silhouette(kmeansObj$cl,dissE^2)
  silinfo<-list(widths=silObj, clus.avg.widths=summary(silObj)$clus.avg.widths, ave.width=summary(silObj)$avg.width)
  return(list(mediods=kmeansObj$centers,clustering=kmeansObj$cluster,call=NA,silinfo=silinfo,objective=NA,diss=dissE,data=x))
}

###Pam
.pamCluster <- function(x,k, ...) { cluster::pam(x=x,k=k,...) }
.pamClassify <- function(x, clusterResult) {
  center<-clusterResult$medoids
  innerProd<-tcrossprod(x,center) #a n x k matrix of inner-products between them
  distMat<-as.matrix(dist(rbind(x,center)))
  distMat<-distMat[1:nrow(x),(nrow(x)+1):ncol(distMat)]
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


