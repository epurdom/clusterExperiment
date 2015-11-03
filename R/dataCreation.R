#' Simulated data for running examples
#'
#' @name simData
#' @aliases simCount
#' @docType data
#' @author Elizabeth Purdom \email{epurdom@@stat.berkeley.edu}
#' @format A data frame with 300 rows and 10 columns
#' @details \code{simData} is simulated normal data of 300 simulated observations with 50 variables, with observations being in one of 3 groups. \code{simCount} is simulated count data of the same dimension. The groups are each of size 100 and are in order in the rows of the data.frames.
#' @keywords data
#' @examples 
#' #code to create data:
#' nvar<-50
#' x<-rbind(matrix(rnorm(100*nvar,mean=5),ncol=nvar),
#' matrix(rnorm(100*nvar,mean=-5),ncol=nvar),
#' matrix(rnorm(100*nvar,mean=0),ncol=nvar))
#' #add in differences in variable means
#' smp<-sample(1:ncol(x),10)
#' x[,smp]<-x[,smp]+10
#' y<-rbind(matrix(rnorm(100*nvar,mean=1),ncol=nvar),
#' matrix(rnorm(100*nvar,mean=-1),ncol=nvar),
#' matrix(rnorm(100*nvar,mean=0),ncol=nvar))
#' y<-y[sample(1:nrow(y)),]+ matrix(rnorm(300*nvar,sd=3),ncol=nvar)
#' simData<-x+y
#' countMean<-exp(simData/2)
#' simCount<-matrix(rpois(n=length(as.vector(countMean)),
#' lambda =as.vector(countMean)+.1),nrow=nrow(countMean),ncol=ncol(countMean))
#' #save(simCount,file="data/simCount.rda")
#' #save(simData,file="data/simData.rda")

NULL