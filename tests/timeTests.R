# code for getting timings of tests:
# Not working for all, just getting "Error: Test failures"
test_out<-testthat::test_dir("testthat",reporter=ListReporter)

results<-sapply(test_out,function(x){
  data.frame(x[c("file","test","real")])
})


# # Moved to longtests
# "105" "test-RSEC.R" "`RSEC` works through whole series of steps" 16.5
# "6" "test-assays.R" "RSEC works independent of assay order" 62.68
# "4" "test-assays.R" "RSEC works wih non default assays" 28.65
# "108" "test-RSEC.R" "`RSEC` works with hdf5" 122.89


# Other possibilities
# "106" "test-RSEC.R" "`RSEC` works with no merging" 7.98
# "107" "test-RSEC.R" "`RSEC` returns clusterMany even when errors later" 2.35

# "76" "test-heatmaps.R" "`plotHeatmap` works with hdf5 objects" 4
# "11" "test-clusterMany.R" "`clusterMany` works with hdf5" 8.53
# "19" "test-clusterSingle.R" "`clusterSingle` works with hdf5Matrix" 10.91

# "73" "test-getBestFeatures.R" "`getBestFeatures` works with weights" 10.11
# "12" "test-clusterMany.R" "`clusterMany` works changing parameters" 8.99
# "54" "test-dendrogram.R" "`makeDendrogram` works with whichCluster" 2.45

# "75" "test-heatmaps.R" "`plotHeatmap` works with matrix objects" 2.37
# "77" "test-heatmaps.R" "`plotHeatmap` works with ClusterExperiment and SummarizedExperiment objects" 14.04
# "78" "test-heatmaps.R" "`plotHeatmap` visualization choices/feature choices all work" 3.87
# "79" "test-heatmaps.R" "`makeBlankData` works" 5.44
# "80" "test-heatmaps.R" "`plotCoClustering` works" 4.53

# "86" "test-mergeClusters.R" "`mergeClusters` works with ClusterExperiment objects" 5.22
# "88" "test-mergeClusters.R" "saving merge info works" 1.66
# "89" "test-mergeClusters.R" "logFC works" 1.92
# "91" "test-mergeClusters.R" "`mergeClusters` works with unassignedSamples" 3.61
#
#
# "93" "test-otherPlots.R" "`plotClusters` works with matrix, ClusterExperiment objects" 2.4
# "94" "test-otherPlots.R" "`plotClusters` rerun above tests with colData included" 2.16
# "95" "test-otherPlots.R" "plotClustersWorkflow" 1.56
# "96" "test-otherPlots.R" "plotting helpers" 4.32
# "98" "test-otherPlots.R" "plotReducedDims works" 2.31
# "100" "test-otherPlots.R" "plotClustersTable works" 22.76
#


# > test_dir(".",reporter=ListReporter)
# Found more than one class "Annotated" in cache; using the first, from namespace 'S4Vectors'
# Also defined by ‘RNeXML’
# Error: Test failures
# > outTests<-test_dir(".",reporter=ListReporter)
# Found more than one class "Annotated" in cache; using the first, from namespace 'S4Vectors'
# Also defined by ‘RNeXML’
# Error: Test failures
# > names(outTests)
# Error: object 'outTests' not found
