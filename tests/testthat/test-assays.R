context("Assays")


test_that("plotting works wih non default assays", {
	## Construct ClusterExperiment with multiple named assays
	suppressMessages(multi_se <- SummarizedExperiment(
	    assays = list(counts = simCount,
	    logcounts = log1p(simCount)))
	    )
	suppressMessages(multi_cc <- ClusterExperiment(multi_se, trueCluster))
  expect_silent(plotFeatureBoxplot(object=multi_cc,
      feature=1,whichAssay=1))
  expect_silent(plotFeatureBoxplot(object=multi_cc,
      feature=1,whichAssay=2))

  suppressMessages(small_cc <- multi_cc[1:10,c(1:3, 290:293)])
  suppressMessages(small_cc <- makeDendrogram(small_cc))
  expect_silent(plotHeatmap(small_cc, whichAssay=1))
  expect_silent(plotHeatmap(small_cc, whichAssay=2))

})


