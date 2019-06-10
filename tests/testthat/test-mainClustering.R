test_that("`mainClustering` options", {

    ########
    #Check mainClustering
    ########
    ###Check pam exactly same:
    expect_silent(x<-mainClustering(mat, 
        inputType="X",
        clusterFunction="pam",
        clusterArgs=list(k=3),
        minSize=1, removeSil=FALSE))
    expect_equal(length(x),ncol(mat))
    x2<-cluster::pam(t(mat),k=3,cluster.only=TRUE)
    expect_equal(x,x2)
    ###Check hierarchicalK exactly same:
    expect_silent(x<-mainClustering(dissMat,
        inputType="diss", 
        clusterFunction="hierarchicalK",
        clusterArgs=list(k=3),
        minSize=1, 
        removeSil=FALSE))
    expect_equal(length(x),ncol(mat))
    x2<-stats::cutree(stats::hclust(as.dist(dissMat)),k=3)
    expect_equal(x,x2)


    #check giving wrong parameters gives warning:
    expect_warning(mainClustering(dissMat, 
        inputType="diss",
        clusterFunction="tight", 
        clusterArgs=list(alpha=0.1),
        minSize=5, removeSil=TRUE),
        "do not match the algorithmType")
    expect_warning(mainClustering(dissMat, inputType="diss",
        clusterFunction="tight", 
        clusterArgs=list(alpha=0.1,k=3),
        minSize=5, 
        removeSil=TRUE),
        "do not match the algorithmType")
    #check missing required parameters gives warning:
    expect_error(mainClustering(dissMat,
        inputType="diss", 
        clusterFunction="tight",
        minSize=5),
        "must supply arguments: alpha")
    expect_error(mainClustering(mat, 
        inputType="X",
        clusterFunction="pam", 
        clusterArgs=list(alpha=0.1),
        minSize=5, removeSil=TRUE),
        "must supply arguments: k")
    #check warning for superfluous arguments passed to clustering function:
    expect_warning(mainClustering(dissMat, 
        inputType="diss",
        clusterFunction="tight", 
        clusterArgs=list(k=3,alpha=0.1),
    	minSize=5),
        "arguments passed via clusterArgs to the clustering function tight are not all applicable")
    expect_warning(mainClustering(mat, 
        inputType="X",
        clusterFunction="pam", 
        clusterArgs=list(k=3,alpha=0.1),
        minSize=5, 
        removeSil=TRUE),
        "arguments passed via clusterArgs to the clustering function pam are not all applicable")
    expect_warning(mainClustering(dissMat, 
        inputType="diss",
        clusterFunction="tight", 
        clusterArgs=list(alpha=0.1, evalClusterMethod="average")),
        "arguments passed via clusterArgs to the clustering function tight are not all applicable")
    expect_warning(mainClustering(dissMat, 
        inputType="diss",
        clusterFunction="hierarchical01", 
        clusterArgs=list(alpha=0.1, minSize.core=4)),
        "arguments passed via clusterArgs to the clustering function hclust are not all applicable")

    #test default 01 distance
    expect_silent(mainClustering(dissMat, inputType="diss",
        clusterFunction="tight", 
        clusterArgs=list(alpha=0.1)))
    #test default K distance
    expect_silent(mainClustering(dissMat, inputType="diss",
        clusterFunction="hierarchicalK", 
        clusterArgs=list(k=3)))

})


test_that("`mainClustering` works with cat", {
    #test default 01 distance
    expect_silent(mainClustering(catMat, inputType="cat",
        clusterFunction="hierarchical01", 
        clusterArgs=list(alpha=0.1)))
    
})

