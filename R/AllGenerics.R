setGeneric(
  name = "clusterAll",
  def = function(x,  subsample=TRUE, sequential=FALSE,
                 clusterFunction=c("tight", "hierarchical", "pam", "kmeans"),
                 clusterDArgs=NULL, subsampleArgs=NULL, seqArgs=NULL) {
    standardGeneric("clusterAll")
  }
)

setGeneric(
  name = "isLog",
  def = function(x) {
    standardGeneric("isLog")
  }
)

setGeneric(
  name = "labels",
  def = function(x) {
    standardGeneric("labels")
  }
)

setGeneric(
  name = "labels<-",
  def = function(object, value) {
    standardGeneric("labels<-")
  }
)
