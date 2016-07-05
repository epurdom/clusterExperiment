#' @rdname plottingFunctions
#' @name plottingFunctions
#' @aliases makeBlankData
#' @param data matrix with samples on columns and features on rows.
#' @param groupsOfFeatures list, with each element of the list containing a
#'   vector of numeric indices.
#' @param nBlankLines the number of blank lines to add in the data matrix to
#'   separate the groups of indices (will govern the amount of white space if
#'   data is then fed to heatmap.)
#' @details \code{makeBlankData} pulls the data corresponding to the row indices
#'   in \code{groupsOfFeatures} adds lines of NA values into data between these
#'   groups. When given to heatmap, will create white space between these groups
#'   of features.
#'
#' @return \code{makeBlankData} returns a list with items
#'  \itemize{
#'  \item{"dataWBlanks"}{ The data with the rows of NAs separating the given
#'  indices.}
#'  \item{"rowNamesWBlanks"}{ A vector of characters giving the rownames for the
#'  data, including blanks for the NA rows. These are not given as rownames to
#'  the returned data because they are not unique. However, they can be given to
#'  the \code{labRow} argument of \code{\link[NMF]{aheatmap}} or
#'  \code{\link{plotHeatmap}}.}
#' }
#'
#' @export
#'
#' @examples
#' data(simData)
#'
#' x <- makeBlankData(simData[,1:10], groupsOfFeatures=list(c(5, 2, 3), c(20,
#' 34, 25)))
makeBlankData <- function(data,groupsOfFeatures,nBlankLines = 1) {
  if (!is.list(groupsOfFeatures))
    stop("groupsOfFeatures must be a list. Will ignore.")

  testIndices <- sapply(groupsOfFeatures,function(x) {
    is.numeric(x) && all(x %in% 1:NROW(data))
  })
  if (!all(testIndices))
    stop("Invalid list of indices in groupsOfFeatures.")
  if(is.null(rownames(data))) row.names(data)<-as.character(1:nrow(data))
  #make list of data of feature groups
  dataList <- lapply(groupsOfFeatures,function(ii){data[ii,,drop=FALSE]})

  #add NA rows between groups
  naData <- matrix(NA,nrow = nBlankLines,ncol = ncol(data))
  dataListMinus <- lapply(dataList[-length(dataList)],function(x) {
    return(rbind(x,naData))
  })
  newData <-
    data.frame(do.call("rbind",c(dataListMinus,dataList[length(dataList)])))

  #make names for this
  rnames <- lapply(dataList,rownames)
  rnamesMinus <-
    lapply(head(rnames,-1),function(x) {
      c(x,rep("",nBlankLines))
    })
  rnames <- unname(c(unlist(rnamesMinus),rnames[[length(rnames)]]))
  #browser()
  #can't set rownames because not unique values!
  #rownames(newData)<-rnames
  return(list(dataWBlanks = newData,rowNamesWBlanks = rnames))

}
