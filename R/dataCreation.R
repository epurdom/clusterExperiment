#' Simulated data for running examples
#'
#' @name simData
#' @aliases simCount 
#' @aliases trueCluster
#' @docType data
#' @author Elizabeth Purdom \email{epurdom@@stat.berkeley.edu}
#' @format Three objects are loaded, two data frame(s) of simulated data each with 300 samples/rows and 50 variables/columns, and a vector of length 300 with the true cluster assignments.
#' @details \code{simData} is simulated normal data of 300 simulated observations with 50 variables, with observations being in one of 3 groups. \code{simCount} is simulated count data of the same dimension. \code{trueCluster} gives the true cluster identifications of the samples. The true clusters are each of size 100 and are in order in the rows of the data.frames. 
#' @keywords data
#' @examples 
#' #code used to create data:
#' nvar<-51 #multiple of 3
#' n<-100
#' x<-rbind(matrix(rnorm(n*nvar,mean=5),ncol=nvar),
#' matrix(rnorm(n*nvar,mean=-5),ncol=nvar),
#' matrix(rnorm(n*nvar,mean=0),ncol=nvar))
#' #make some of them flipped effects (better for testing if both sig under/over expressed variables)
#' geneGroup<-sample(rep(1:3,each=floor(nvar/3)))
#' gpIndex<-list(1:n,(n+1):(n*2),(2*n+1):(n*3))
#' x[,geneGroup==1]<-x[unlist(gpIndex[c(3,1,2)]),geneGroup==1]
#' x[,geneGroup==2]<-x[unlist(gpIndex[c(2,3,1)]),geneGroup==2]
#' 
#' #add in differences in variable means
#' smp<-sample(1:ncol(x),10)
#' x[,smp]<-x[,smp]+10
#' y<-rbind(matrix(rnorm(n*nvar,mean=1),ncol=nvar),
#' matrix(rnorm(n*nvar,mean=-1),ncol=nvar),
#' matrix(rnorm(n*nvar,mean=0),ncol=nvar))
#' y<-y[sample(1:nrow(y)),]+ matrix(rnorm(3*n*nvar,sd=3),ncol=nvar)
#' simData<-x+y
#' #add pure noise variables
#' simData<-cbind(simData,matrix(rnorm(3*n*nvar,mean=10),ncol=nvar),matrix(rnorm(3*n*nvar,mean=5),ncol=nvar))
#' countMean<-exp(simData/2)
#' simCount<-matrix(rpois(n=length(as.vector(countMean)),
#' lambda =as.vector(countMean)+.1),nrow=nrow(countMean),ncol=ncol(countMean))
#' trueCluster<-rep(c(1:3),each=n)
#' #save(list=c("simCount","simData","trueCluster"),file="../data/simData.rda")
NULL
