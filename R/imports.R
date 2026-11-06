# Package-wide imports of base-R and Biobase functions used across many files.
# Before NMF stopped being imported (see plotHeatmap.R), these resolved only
# because NMF's .onLoad attached Biobase (and through it BiocGenerics, which
# attaches stats and graphics). They must be imported explicitly.
#' @importFrom graphics abline axis par plot.new plot.window rect strwidth
#'   text title
#' @importFrom stats as.dist cor mad median model.matrix na.omit nobs p.adjust
#'   rnorm runif var
#' @importFrom utils packageVersion
#' @importFrom Biobase exprs
NULL
