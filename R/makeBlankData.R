# markers <- read.table("../data/allen_markers.txt", as.is=TRUE)
# allen_genes <- intersect(markers[,1], rownames(norm))
# allen_clusters <- markers[,2]
# names(allen_clusters) <- markers[,1]
# allen_clusters <- as.factor(allen_clusters[allen_genes])
# 
# only_4 <- which(final_merged %in% c(11, 3, 8, 12, 9))
# names(final_merged) <- colnames(norm)
# cl <- droplevels(final_merged[only_4])
# cl <- factor(cl, levels=c(11, 3, 8, 12, 9))
# cl <- sort(cl)
# keep <- names(cl)
# 
# genes1 <- c("Etv1", "Fam84b", "Tcerg1l", "Crym")
# genes2 <- c("Cdh13", "Syt17", "Cpne7", "Nnat")
# genes3 <- c("Endou", "Whrn", "Ddit4l", "Aldh1l1", "Plb1")
# genes4 <- c("Deptor", "Myl4", "Tmem91", "Il1rapl2", "Hsd11b1")
# genes5 <- c("Mc4r", "Ctxn3", "Scml2")
# 
# geneord <- c(genes1, genes2, genes3, genes4)
# groupLength <- c(length(genes1), length(genes2), length(genes3), length(genes4))
# sep <- cumsum(groupLength)
# 
# #function to add blank lines of NA in data of arbitrary size
# makeBlankData<-function(data,sep,nadd=1){
#   naData<-matrix(NA,nrow=nadd,ncol=ncol(data))
#   colnames(naData)<-colnames(data)
#   start<-c(1,sep+1)
#   end<-c(sep,nrow(data))
#   len<-end-start+1
#   grFac<-rep(1:length(len),times=len)
#   dataList<-by(data,grFac,function(x){x})
#   dataListMinus<-lapply(dataList[-length(dataList)],function(x){
#     return(rbind(x,naData))
#   })
#   rnames<-lapply(dataList,rownames)
#   rnamesMinus<-lapply(head(rnames,-1),function(x){c(x,rep("",nadd))})
#   rnames<-c(unlist(rnamesMinus),rnames[[length(rnames)]])
#   return(list(data=do.call("rbind",c(dataListMinus,dataList[length(dataList)])),rownames=rnames))
# }
# 
# temp<-makeBlankData(norm[geneord,keep],sep=head(sep,-1),nadd=2)
# dat<-t(temp$data)
# rnames<-temp$rownames
# 
# out <- dualHeatmap(clusterVec=cl, heatData=t(norm[,keep]), dual=FALSE, clusterSamples=FALSE, clusterVars=FALSE, annCol=data.frame(Cluster=cl), annColors=list(Cluster=colMerged[c(7, 2, 5, 8)]), main="", whVars = geneord, annLegend=FALSE)
# 
# dualHeatmap(clusterVec=cl, heatData=dat, dual=FALSE, clusterSamples=FALSE, clusterVars=FALSE, annCol=data.frame(Cluster=cl), annColors=list(Cluster=colMerged[c(7, 2, 5, 8, 6)]), main="Allen Institute L5 markers", labRow=rnames, annLegend=FALSE, breaks=out$breaks)
