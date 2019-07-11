context("RsecFluidigm")

test_that("saved rsecFluidigm is still valid object", {
	data(rsecFluidigm)
	expect_true(validObject(rsecFluidigm))
    
    #checks that current code will still recreate same basic rsecFluidigm
    # takes a bit of time, so skip on bioc, just want on travis so get alert
    # that need to deal with this.
    skip_on_bioc()
    rsecFluidigmNew<-makeRsecFluidigmObject(se)
    checkRsecFluidigmObject(rsecFluidigmNew)
    
})