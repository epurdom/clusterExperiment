context("RsecFluidigm")

test_that("saved rsecFluidigm is still valid object", {
	data(rsecFluidigm)
	expect_true(validObject(rsecFluidigm))
    
    #checks that current code will still recreate same basic rsecFluidigm
    # takes a bit of time, so skip on bioc, just want on travis so get alert
    # that need to deal with this.
    skip_on_bioc()
    data(fluidigmData)
    data(fluidigmColData)
    expect_silent(se<-SummarizedExperiment(assays=fluidigmData,
        colData=fluidigmColData))
    expect_message(rsecFluidigmNew<-makeRsecFluidigmObject(se),
        "Will calculate all the methods requested")
    expect_silent(checkRsecFluidigmObject(rsecFluidigmNew))
    data(rsecFluidigm)
    expect_equal(clusterMatrix(rsecFluidigmNew),
        clusterMatrix(rsecFluidigm))
    
    
})
