library(clusterCells)

mat <- matrix(data=rnorm(200), ncol=10)
mat[1,1]<- -1 #force a negative value
labels <- gl(5, 2)
se <- SummarizedExperiment(mat)

cc <- clusterCells(mat, as.numeric(labels), transformation = function(x){x})
cc2 <- clusterCells(se, as.numeric(labels), transformation = function(x){x})

test_that("`clusterCells` constructor works with matrix and
          SummarizedExperiments", {
            expect_equal(cc, cc2)
            expect_error(clusterCells(mat), "missing")
            expect_error(clusterCells(mat,as.numeric(labels),transformation=log),info="Error checking transFun")
            expect_error(clusterCells(mat, as.numeric(labels)), "missing")
            expect_error(clusterCells(mat, labels[1:2], function(x){x}),
                         "must be a vector of length equal to the number of samples")
            expect_error(clusterCells(as.data.frame(mat), labels, function(x){x}),
                         "must be a matrix or SummarizedExperiment object")
            expect_error(clusterCells(mat, as.character(labels), function(x){x}),
                         "must be a numeric matrix")

            expect_warning(clusterCells(mat, labels, function(x){x}), "was coerced to numeric")

            expect_is(cc, "ClusterCells")
            expect_is(cc, "SummarizedExperiment")
          })

test_that("adding clusters, setting primary labels and remove unclustered cells
          work as promised", {
            c1 <- addClusters(cc, rep(c(-1, 1), each=5))
            expect_equal(NCOL(allClusters(cc)), 1)
            expect_equal(NCOL(allClusters(c1)), 2)
            expect_equal(primaryCluster(c1), primaryCluster(cc))
            expect_is(transformation(cc),"function")
            primaryClusterIndex(c1) <- 2
            expect_false(all(primaryCluster(c1)==primaryCluster(cc)))

            expect_equal(cc, removeUnclustered(cc))

            c2 <- removeUnclustered(c1)
            expect_equal(NCOL(c2), 5)
          })
