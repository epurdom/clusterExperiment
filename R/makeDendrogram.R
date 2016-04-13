#' @title Make hierarchy of set of clusters
#'
#' @description Makes a dendrogram of a set of clusters based on hclust on the mediods of the cluster.
#'
#'
#' @param x data to define the mediods from. Samples on rows, features on columns
#' @param cluster A numeric vector with cluster assignments. If x is a ClusterExperiment class, cluster is automatically the primaryCluster(x). ``-1'' indicates the sample was not assigned to a cluster.
#' @param leaves character value to indicate what should be the leaves. See details.
#' @param unassigned how to handle unassigned samples("-1") ; only relevant if leaves="samples" See details. 
#' @param ... arguments passed to hclust
#' @details The dendrogram is made based on the cluster assignments given in cluster. If leaves="clusters", the dendrogram is created by applying \code{\link{hclust}} to the mediods of each cluster. If \code{leaves="samples"} then the clusters are again clustered, but now the samples are also part of the resulting dendrogram. This is done by giving each sample the value of the mediod of its cluster. 
#' 
#' @details The argument \code{unassignedSamples} governs what is done with unassigned samples (defined by a -1 cluster value). If unassigned=="cluster", then the dendrogram is created by hclust of the expanded mediod data plus the original unclustered observations. If \code{unassignedSamples} is "outgroup", then all unassigned samples are put as an outgroup in 
#'
#'
#' @return A dendrogram for the clusters. If full=TRUE, the number of leaves is equal to the number of samples. If full=FALSE, the number of leaves is equal to the number of clusters.

#' @examples
#' data(simData)
#' #create a clustering, for 8 clusters (truth was 3)
#' cl<-clusterSingle(simData,clusterFunction="pam",subsample=FALSE,
#' sequential=FALSE, clusterDArgs=list(k=8))$cl
#' #create dendrogram of clusters:
#' hcl<-clusterHclust(dat=simData,cl,full=FALSE)
#' plot(hcl)
#'
#' #create dendrogram for plotting with data in heatmap:
#' hclData<-clusterHclust(dat=simData,cl,full=TRUE)
#' plotHeatmap(cl,heatData=simCount,clusterData=hclData,colorScale=seqPal5,
#'	annCol=data.frame(PAM8=cl,Truth=trueCluster))
#' @rdname makeDendrogram
setMethod(
  f = "makeDendrogram",
  signature = "ClusterExperiment",
  definition = function(x,...)
  {
    makeDendrogram(x=assay(x),cl=primaryCluster(x),...)
  })

#' @rdname makeDendrogram
setMethod(
  f = "makeDendrogram",
  signature = "matrix",
  definition = function(x,cluster,
                   leaves=c("samples","clusters"),
                   unassignedSamples=c("outgroup","cluster","remove"),...){
    leaves<-match.arg(leaves)
    unassigned<-match.arg(unassignedSamples)
    dat<-x
    cl<-cluster
    if(leaves=="samples") full<-TRUE else full<-FALSE
    if(ncol(dat)!=length(cl)) stop("cl must be the same length as the number of columns of dat")
    if(is.null(colnames(dat))) colnames(dat)<-as.character(1:ncol(dat))
    if(is.factor(cl)){warning("cl is a factor. Converting to numeric, which may not result in valid conversion")
      cl<-as.numeric(as.character(cl))
    }
    dat<-t(dat) #make like was in old code
    
    #############
    whRm<- which(cl>=0) #remove -1, -2
    if(length(whRm)==0){
        warning("all samples have clusterIds <0. Will just use all clustersIds")
        whRm<-1:length(cl)
    }
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
  
})

#' @rdname plotDendrogram
#' @inheritParams makeDendrogram
setMethod(
    f = "plotDendrogram",
    signature = "ClusterExperiment",
    definition = function(x,leaves=c("samples","clusters"),...)
    {
        leaves<-match.arg(leaves)
        if(leaves=="samples") plot(dendrogram(x),...)
        else{plot(makeDendrogram(x,leaves="clusters"),...)}
    })
