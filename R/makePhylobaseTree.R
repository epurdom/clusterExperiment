####
#Convert to object used by phylobase so can navigate easily -- might should make generic function...
.makePhylobaseTree<-function(x,type){
	type<-match.arg(type,c("hclust","dendro"))
	if(type=="hclust"){
		#first into phylo from ape package
		tempPhylo<-try(ape::as.phylo(x),FALSE)
		if(inherits(tempPhylo, "try-error")) stop("the hclust object cannot be converted to a phylo class with the methods of the 'ape' package.")
	}
	if(type=="dendro"){
		tempPhylo<-try(dendextend::as.phylo.dendrogram(x),FALSE)
		if(inherits(tempPhylo, "try-error")) stop("the dendrogram object cannot be converted to a phylo class with the methods of 'dendextend' package. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
	}
	phylo4Obj<-try(as(tempPhylo,"phylo4"),FALSE) 
	if(inherits(phylo4Obj, "try-error")) stop("the internally created phylo object cannot be converted to a phylo4 class. Check that you gave simple hierarchy of clusters, and not one with fake data per sample")
	phylobase::nodeLabels(phylo4Obj)<-paste("Node",1:phylobase::nNodes(phylo4Obj),sep="")
	return(phylo4Obj)
}