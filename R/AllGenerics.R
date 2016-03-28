setGeneric(
  name = "clusterAll",
  def = function(x,  subsample=TRUE, sequential=FALSE,
                 clusterFunction=c("tight", "hierarchical", "pam", "kmeans"),
                 clusterDArgs=NULL, subsampleArgs=NULL, seqArgs=NULL) {
    standardGeneric("clusterAll")
  }
)
