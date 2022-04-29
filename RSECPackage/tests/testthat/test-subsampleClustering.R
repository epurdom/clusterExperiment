context("subsampleClustering")

test_that("subsampling", {

    set.seed(4897)
    expect_silent(subAll <- subsampleClustering(inputMatrix=mat, 
        inputType="X",clusterFunction="pam",
        clusterArgs=list(k=3), classifyMethod="All",
        resamp.num = 100,
        samp.p = 0.7,ncores=1))




	#subsample clusterings won't have identification to all samples...
    set.seed(4897)
    expect_silent(subInSample <- subsampleClustering(inputMatrix=mat,
            inputType="X",clusterArgs=list(k=3),
            clusterFunction="pam",
            classifyMethod=c("InSample"),
            resamp.num = 100,
            samp.p = 0.7,ncores=1))

	
    ## Windows does not support mclapply
    skip_on_os("windows")
	set.seed(4897)
	expect_silent(subAllParal <- 
	    subsampleClustering(inputMatrix=mat,
            inputType="X",clusterFunction="pam",
	        clusterArgs=list(k=3), classifyMethod="All",
	        resamp.num = 100,samp.p = 0.7,ncores=2)
	)

	expect_identical(subAllParal,subAll)

	

})

