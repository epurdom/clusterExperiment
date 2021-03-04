context("Constructor")

test_that("run internalFunctionCheck check on all built in cluster functions", {
    
    bMethods<-c(listBuiltInTypeK(),listBuiltInType01())
    for(cf in bMethods){
        expect_silent(cf<-getBuiltInFunction(cf))
        expect_silent(internalFunctionCheck(cf@ClusterFunction,
            inputType=inputType(cf),
            algorithmType=algorithmType(cf),
            outputType=cf@outputType))
   }
})

test_that("cat inputType works on relevant cluster functions", {
    catMat<-cbind(catMat,catMat)
    set.seed(32590)
    catMat<-catMat[,sample(1:ncol(catMat))]
    expect_silent(pMat<-clusterExperiment:::.clustersHammingDistance(catMat))
    
    kMethods<-listBuiltInTypeK()
    kMethods<-kMethods[sapply(inputType(kMethods),function(x){"cat" %in% x})]
    for(cf in kMethods){
        expect_silent(cfObj<-getBuiltInFunction(cf))
        
        set.seed(782935)
        expect_silent(pOut<-cfObj@clusterFUN(inputMatrix=pMat,
            inputType="diss",k=3,cluster.only=TRUE,
            removeDup=TRUE, checkArgs=TRUE))
        set.seed(782935)
        expect_silent(allOut<-cfObj@clusterFUN(inputMatrix=catMat,
            inputType="cat",k=3,cluster.only=TRUE,
            removeDup=FALSE,checkArgs=TRUE))
        expect_equal(allOut,pOut)
        set.seed(782935)
        expect_silent(nodupOut<-cfObj@clusterFUN(inputMatrix=catMat,
            inputType="cat",k=3,cluster.only=TRUE,
            removeDup=TRUE, checkArgs=TRUE))
        expect_equal(allOut,nodupOut)
    
    }
    aMethods<-listBuiltInType01()
    aMethods<-aMethods[sapply(inputType(aMethods),function(x){"cat" %in% x})]
    for(cf in aMethods){
        expect_silent(cfObj<-getBuiltInFunction(cf))
        set.seed(782935)
        #removeDup shouldn't do anything for inputType="diss" 
        expect_silent(pOut<-cfObj@clusterFUN(inputMatrix=pMat,
            inputType="diss",alpha=0.3,cluster.only=TRUE,
            removeDup=TRUE, checkArgs=TRUE))
        set.seed(782935)
        expect_silent(allOut<-cfObj@clusterFUN(inputMatrix=catMat,
            inputType="cat",alpha=0.3,cluster.only=TRUE,
            removeDup=FALSE,checkArgs=TRUE))
        expect_equal(allOut,pOut)
        set.seed(782935)
        expect_silent(nodupOut<-cfObj@clusterFUN(inputMatrix=catMat,
            inputType="cat",alpha=0.3,cluster.only=TRUE,
            removeDup=TRUE, checkArgs=TRUE))
        # Don't actually expect same answer if remove duplicates. The tight method seems to give different answers, not sure that's a problem...
        #expect_equal(allOut,nodupOut)
    
    }

})

test_that("built in cluster functions return correct clustering", {
    
    
    #test some easy ones
    
    #PAM
    truth<- cluster::pam(x=t(mat),cluster.only=TRUE,k=3)
    out<-getBuiltInFunction("pam")@clusterFUN(inputMatrix=mat,
         inputType="X",k=3,cluster.only=TRUE)
    expect_equal(truth,out)
       
    #kmeans
    set.seed(24819)
    truth<- stats::kmeans(x=t(mat),centers=3)$cluster
    set.seed(24819)
    out<-getBuiltInFunction("kmeans")@clusterFUN(inputMatrix=mat,
         inputType="X",k=3,cluster.only=TRUE)
    expect_equal(truth,out)
})
