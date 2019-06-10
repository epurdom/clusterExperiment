context("Constructor")

test_that("run ClusterFunction check on all built in functions", {
    
    kMethods<-listBuiltInTypeK()
    for(cf in kMethods){
        expect_silent(internalFunctionCheck(cf@ClusterFunction,
            inputType=inputType(cf),
            algorithmType=algorithmType(cf),
            outputType=cf@outputType))
   }
   aMethods<-listBuiltInType01()
   for(cf in aMethods){
       expect_silent(internalFunctionCheck(cf@ClusterFunction,
           inputType=inputType(cf),
           algorithmType=algorithmType(cf),
           outputType=cf@outputType))
  } 
})

test_that("built in cluster functions work", {
    
    
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
