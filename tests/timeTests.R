# code for getting timings of tests:
# Not working for all, just getting "Error: Test failures"
test_out<-testthat::test_dir("testthat",reporter=ListReporter)

results<-sapply(test_out,function(x){
  data.frame(x[c("file","test","real")])
})


# # Consider moving to longtests:
# "105" "test-RSEC.R" "`RSEC` works through whole series of steps" 16.5
# "11" "test-clusterMany.R" "`clusterMany` works with hdf5" 8.53
# "6" "test-assays.R" "RSEC works independent of assay order" 62.68
# "4" "test-assays.R" "RSEC works wih non default assays" 28.65