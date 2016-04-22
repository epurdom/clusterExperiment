context("Dendrogram")
source("create_objects.R")

mat <- matrix(data=rnorm(200), ncol=10)
mat[1,1]<- -1 #force a negative value
labels <- as.character(gl(5, 2))
labels[c(1:2)]<- c("-1","-2") #make sure some not assigned
labels<-factor(labels)
labMat<-cbind(as.numeric(as.character(labels)),as.numeric(as.character(labels)))
se <- SummarizedExperiment(mat)
cc2 <- clusterExperiment(se, labMat, transformation = function(x){x})
test_that("`makeDendrogram` works with matrix, ClusterExperiment objects", {
    #test matrix version
    makeDendrogram(mat,primaryCluster(cc2))
    #test CE version
    makeDendrogram(cc2)
    #test CE version
    makeDendrogram(cc2,unassigned="cluster")
    
    #test matrix version
    expect_equal(nobs(makeDendrogram(mat, primaryCluster(cc2), unassigned="remove")$samples),length(labels)-2)
    #add test for CE version
})
