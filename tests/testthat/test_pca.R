context("PCA")


test_that("Fast PCA gives the same results as PCA", {

  ## k = NCOL(mat) (should use regular svd)
  pca_res <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  expect_warning(pca_res2 <- clusterExperiment:::.pcaDimRed(mat, md=NCOL(mat), isPct = FALSE, rowvars = rowVars(mat)),"all singular values are requested")
  expect_equivalent(pca_res$x, pca_res2)

  ## k < NCOL(mat) -- note that the signed of some components may be flipped and irlba is an approximation
  pca_res <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  pca_res2 <- clusterExperiment:::.pcaDimRed(mat, md=10, isPct = FALSE, rowvars = rowVars(mat))
  expect_equivalent(abs(pca_res$x[,1:10]), abs(pca_res2))

  pca_res <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  pca_res2 <- clusterExperiment:::.pcaDimRed(mat, md=2, isPct = FALSE, rowvars = rowVars(mat))
  expect_equivalent(abs(pca_res$x[,1:2]), abs(pca_res2), tolerance = 1e-5)

})
