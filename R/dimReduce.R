#' General wrapper to perform dimensionality reduction
#'
#' Given a \code{\link{ClusterExperiment}} object, it returns a copy of the same
#' object with the additional dimensionality reduced data in the
#' \code{dimReduce} slot.
#'
#' Note that with this function one can perform several dimensionality reduction
#' techniques on the data at the same time. The results will be stored as a
#' \code{SimpleList} in the \code{reducedDim} slot.
#'
#' For now, only PCA is implemented.
#'
#' @param x a \code{SingleCellExperiment} object.
#' @param methods a vector of character strings that specify which methods to
#'   use for dimensionality reduction.
#' @param k the number of dimensions.
#'
#' @export
#' @rdname dimReduction
setMethod(
  f = "dimReduction",
  signature = "ClusterExperiment",
  definition = function(x, methods = "PCA", k) {

    method <- match.arg(method)

    if("PCA" %in% methods) {
      x <- pca(x, k)
    }

    return(x)

  }
)

#' Wrapper to perform PCA on a SingleCellExperiment
#'
#' Given a \code{\link{ClusterExperiment}} object, it
#' returns a copy of the same object with the principal components added in the
#' \code{dimReduce} slot with the name "PCA".
#'
#' This method uses the \code{\link[RSpectra]{svds}} function which is faster
#' than \code{svd} when \code{k} is small.
#'
#' @param x a \code{SingleCellExperiment} object.
#' @param k the number of principal components.
#' @param center a logical value indicating whether the variables should be
#'   shifted to be zero centered.
#' @param scale a logical value indicating whether the variables should be
#'   scaled to have unit variance.
#'
#' @importFrom RSpectra svds
#'
#' @export
#' @rdname pca
setMethod(
  f = "pca",
  signature = "ClusterExperiment",
  definition = function(x, k, center = TRUE, scale = FALSE) {

    t_x <- t(transform(x))
    svd_raw <- svds(scale(t_x, center=center, scale=scale), k=k, nu=k, nv=0)
    pc_raw <- svd_raw$u %*% diag(svd_raw$d, nrow = length(svd_raw$d))
    colnames(pc_raw) <- paste0("PC", seq_len(k))
    rownames(pc_raw) <- colnames(x)
    reducedDim(x, "PCA") <- pc_raw
    return(x)
  }
)

## this function should compute the most variable genes and store it as
## (internal?) rowData this is useful as a way to store the data without
## increasing the size of the object but the user should have the illusion that
## the reduced dataset exists and there will be another function
## (getMostVariableFeatures?) used by the user to access the data this function
## is used in clusterSingle et al., but should probably not be exported (?)
# setMethod(
#   f = "computeMostVariableFeatures",
#   signature = "ClusterExperiment",
#   definition = function(x, method = c("var", "mad", "cv")) {
#
#   }
# )
