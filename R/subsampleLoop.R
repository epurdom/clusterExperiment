#---------
#for a sample index ii, determine what other indices in the same sample with it
#ii an index of a sample
#clustVec the clusterIds as returned by above in DList
#clusterLengths the length of each cluster as returned by above in DList
otherIds<-function(idx,clustVec,clustLeng){
    m<-which(clustVec==idx)
    if(length(m)>1) stop("ids clustered in more than one cluster")
    if(length(m)==0) return(NA) #sample not ever clustered
    if(length(m)==1){
        ends<-cumsum(clustLeng)
        begins<-cumsum(c(1,head(clustLeng,-1)))
        whCluster<-which(m<=ends & m>=begins)
        if(length(whCluster)>1 | length(whCluster)==0) stop("error in coding: finding range of clusterids")
        return(clustVec[seq(begins[whCluster],ends[whCluster],by=1)])
    }
}
#Test: otherIds(5,DList[[1]][[1]],DList[[1]][[2]])

#---------
#For ii an index of a sample, calculates the proportion of times joint with every other sample with index < ii.
#clusterList is the results of subsampling, i.e. list with indices of clusters adjacent
#returns vector of length N-ii with the proportions
#NA if either never subsampled together or never clustered together.
#---------
searchForPairs<-function(ii,clusterList,N){
    #get list of those indices sample ii was sampled
	whHave<-which(sapply(clusterList,function(ll){ii%in%ll$clusterIds}))
	#calculate number of times sampled with (denominator)
    sampledWithTab<-table(unlist(sapply(clusterList[whHave],.subset2,"clusterIds")))
    #get those indices clustered with and tabulate
	clusterWith<-lapply(clusterList[whHave],function(ll){
        otherIds(idx=ii,clustVec=ll$clusterIds,clustLeng=ll$clusterLengths)
    })
	# if(doGC){#just to reduce memory since will be parallelized
	# 	rm(clusterList)
	# 	gc()
	# }
	clusterWithTab<-table(unlist(clusterWith))
    jointNames<-as.character(1:N)
	whLower<-which(as.integer(as.numeric(jointNames))<ii)
	return(as.integer(clusterWithTab[jointNames][whLower])/as.integer(sampledWithTab[jointNames][whLower]))

	#thoughts about alternative code if not need save NxN matrix...
    #jointNames<-names(sampledWithTab) #if manage to not save NxN matrix, could use this to return only those that actually present
    #out<-out[!is.na(out[,"together"]),,drop=FALSE] #if manage to not save NxN matrix, could use this to return only those that actually present; but then need to not return proportions, but something else.
}
