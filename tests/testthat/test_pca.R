context("PCA")
source("create_objects.R")

test_that("Fast PCA gives the same results as PCA", {

  ## k = NCOL(mat) (should use regular svd)
  pca_res <- prcomp(mat)
  expect_warning(pca_res2 <- clusterExperiment:::.pca(mat, k=NCOL(mat)),"all singular values are requested")
  expect_equivalent(pca_res$x, pca_res2)

  pca_res <- prcomp(mat, center=FALSE)
  expect_warning(pca_res2 <- clusterExperiment:::.pca(mat, k=NCOL(mat), center=FALSE),"all singular values are requested")
  expect_equivalent(pca_res$x, pca_res2)

  pca_res <- prcomp(mat, scale=TRUE)
  expect_warning(pca_res2 <- clusterExperiment:::.pca(mat, k=NCOL(mat), scale=TRUE),"all singular values are requested")
  expect_equivalent(pca_res$x, pca_res2)

  pca_res <- prcomp(mat, scale=TRUE, center=FALSE)
  expect_warning(pca_res2 <- clusterExperiment:::.pca(mat, k=NCOL(mat), scale=TRUE, center=FALSE),"all singular values are requested")
  expect_equivalent(pca_res$x, pca_res2)

  ## k < NCOL(mat) -- note that the signed of some components may be flipped
  pca_res <- prcomp(mat)
  pca_res2 <- clusterExperiment:::.pca(mat, k=10)
  expect_equivalent(abs(pca_res$x[,1:10]), abs(pca_res2))

  pca_res <- prcomp(mat, center=FALSE)
  pca_res2 <- clusterExperiment:::.pca(mat, k=10, center=FALSE)
  expect_equivalent(abs(pca_res$x[,1:10]), abs(pca_res2))

  pca_res <- prcomp(mat, center=FALSE, scale=TRUE)
  pca_res2 <- clusterExperiment:::.pca(mat, k=10, center=FALSE, scale=TRUE)
  expect_equivalent(abs(pca_res$x[,1:10]), abs(pca_res2))

  pca_res <- prcomp(mat, center=TRUE, scale=TRUE)
  pca_res2 <- clusterExperiment:::.pca(mat, k=10, center=TRUE, scale=TRUE)
  expect_equivalent(abs(pca_res$x[,1:10]), abs(pca_res2))

})
