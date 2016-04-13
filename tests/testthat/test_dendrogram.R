mat <- matrix(data=rnorm(200), ncol=10)
mat[1,1]<- -1 #force a negative value
labels <- as.character(gl(5, 2))
labels[c(1:2)]<- c("-1","-2") #make sure some not assigned
labels<-factor(labels)
labMat<-cbind(labels,labels)
se <- SummarizedExperiment(mat)
cc2 <- clusterExperiment(se, as.numeric(labels), transformation = function(x){x})
test_that("`makeDendrogram` works with matrix, ClusterExperiment objects", {
    #test matrix version
    makeDendrogram(mat,primaryCluster(cc2),leaves="samples")
    #test CE version
    makeDendrogram(cc2)
    #test CE version
    makeDendrogram(cc2,leaves="clusters")
    #test CE version
    makeDendrogram(cc2,leaves="samples",unassigned="cluster")
    #test CE version
    expect_equal(makeDendrogram(cc2,leaves="samples",unassigned="remove"),length(labels)-2)
    #test assignment
    dendrogram(cc2)<-makeDendrogram(cc2)
    expect_error(dendrogram(cc2)<-makeDendrogram(cc2,leaves="clusters"))
    })
