#' @title Create contrasts for testing DE of a cluster
#' @docType methods
#' @description Uses clustering to create different types of contrasts to be
#'   tested that can then be fed into DE testing programs.
#' @aliases clusterContrasts
#' @param cluster Either a vector giving contrasts assignments or a
#'   ClusterExperiment object
#' @param contrastType What type of contrast to create.
#'   `Dendro' traverses the given dendrogram and does contrasts of the samples
#'   in each side,  `Pairs' does pair-wise contrasts based on the pairs given in
#'   pairMat (if pairMat=NULL, does all pairwise), and `OneAgainstAll' compares
#'   each cluster to the average of all others.
#' @param dendro The dendrogram to traverse if contrastType="Dendro". Note that this
#'   should be the dendrogram of the clusters, not of the individual samples.
#' @param pairMat matrix giving the pairs of clusters for which to do pair-wise
#'   contrasts (must match to elements of cl). If NULL, will do all pairwise of
#'   the clusters in \code{cluster} (excluding "-1" categories). Each row is a pair
#'   to be compared and must match the names of the clusters in the vector
#'   \code{cluster}.
#' @param removeNegative logical, whether to remove negative valued clusters 
#'   from the design matrix. Appropriate to pick TRUE (default) if design will
#'   be input into linear model on samples that excludes -1.
#' @param outputType character string. Gives format for the resulting contrast
#'   matrix. Currently the only option is the format appropriate for
#'   \code{\link{limma}} package, but we anticipate adding more.
#' @param ... arguments that are passed to from the \code{ClusterExperiment}
#'   version to the most basic numeric version.
#' @details The input vector must be numeric clusters, but the external commands
#'   that make the contrast matrix (e.g. \code{\link{makeContrasts}}) require
#'   syntatically valid R names. For this reason, the names of the levels will
#'   be "X1" instead of "1". And negative values (if removeNegative=FALSE) will
#'   be "X.1","X.2", etc.
#' @return List with components:
#'\itemize{
#' \item{\code{contrastMatrix}}{ Contrast matrix, the form of which depends on \code{outputType}.
#' If \code{outputType=="limma"}, the result of running
#'   \code{\link{makeContrasts}}: a matrix with number of columns equal
#'   to the number of contrasts, and rows equal to the number of levels of the
#'   factor that will be fit in a linear model.}
#'   \item{\code{contrastNames}}{A vector of names for each of the contrasts. NULL if no such additional names.}
#'   }
#' @author Elizabeth Purdom
#' @examples 
#' data(simData)
#' cl <- clusterMany(simData,nPCADims=c(5,10,50),  dimReduce="PCA",
#' clusterFunction="pam", ks=2:4, findBestK=c(FALSE), removeSil=TRUE,
#' subsample=FALSE)
#' #Pairs:
#' clusterContrasts(cl,contrastType="Pairs")
#' #Dendrogram
#' cl<-makeDendrogram(cl)
#' clusterContrasts(cl,contrastType="Pairs")
#' @export
#' @rdname clusterContrasts
setMethod(f = "clusterContrasts",
          signature = "ClusterExperiment",
          definition = function(cluster,contrastType,...){
            if(contrastType=="Dendro"){
              if(is.null(cluster@dendro_clusters)) stop("Must run makeDendrogram before calling clusterContrasts if want contrastType='Dendro'")
              else dendro<-cluster@dendro_clusters
            }
            else dendro<-NULL
            clusterContrasts(primaryCluster(cluster),contrastType=contrastType,dendro=dendro,...)
          })
#' @rdname clusterContrasts
#' @export
#' @importFrom limma makeContrasts
#' @importFrom MAST Hypothesis
setMethod(f = "clusterContrasts",
          signature = "vector",
          definition = function(cluster,contrastType=c("Dendro", "Pairs", "OneAgainstAll"),
                                dendro=NULL, pairMat=NULL,outputType="limma",removeNegative=TRUE){
            cluster<-.convertToNum(cluster)
           if(removeNegative) cl<-cluster[cluster>0] else cl<-cluster
            cl<-factor(cl)
            contrastType<-match.arg(contrastType)
            if(contrastType=="Dendro"){
              if(is.null(dendro)) stop("must provide dendrogram if contrastType='Dendro'")
              ####
              #Convert to object used by phylobase so can navigate easily -- might should make generic function...
              phylo4Obj<-.makePhylobaseTree(dendro,type="dendro")
              clChar<-as.character(cl)
              allTipNames<-phylobase::labels(phylo4Obj)[phylobase::getNode(phylo4Obj,  type=c("tip"))]
              if(any(sort(allTipNames)!=sort(unique(clChar)))) stop("tip names of dendro don't match cluster vector values")
              
              #each internal node (including root) construct contrast between them.
              #(before just tested differences between them)
              allInternal<-phylobase::getNode(phylo4Obj,  type=c("internal"))
              .makeNodeContrast<-function(nodeId){
                children<-phylobase::descendants(phylo4Obj,nodeId,"children") #get immediate descendants
                if(length(children)!=2) stop("More than 2 children for internal node; does not make sense with code")
                #find tips of each child:
                contrAvePerChild<-sapply(children,function(x){
                  tips<-phylobase::descendants(phylo4Obj,x,"tip")
                  tipNames<-phylobase::labels(phylo4Obj)[tips]
                  #should make this code use make.names instead of pasting X...
                  if(length(tipNames)>1) return(paste("(",paste(make.names(tipNames),collapse="+"),")/",length(tips),sep=""))
                  else return(make.names(tipNames))
                })
                return(paste(contrAvePerChild,collapse="-"))
              }
              contrastNames<-sapply(allInternal,.makeNodeContrast)
              
            }
            if(contrastType=="OneAgainstAll"){
              levs<-levels(cl)
              contrastNames<-sapply(levs,function(x){
                one<-make.names(x)
                all<-make.names(levs[-which(levs==x)])
                all<-paste("(",paste(all,collapse="+"),")/",length(all),sep="")
                contr<-paste(one,"-",all,sep="")
                return(contr)
              })
              names(contrastNames)<-levs
            }
            if(contrastType=="Pairs"){
              if(is.null(pairMat)){ #make pair Mat of all pairwise
                levs<-levels(cl)
                pairMat<-t(apply(expand.grid(levs,levs),1,sort))
                pairMat<-unique(pairMat)
                pairMat<-pairMat[-which(pairMat[,1]==pairMat[,2]),,drop=FALSE]
              }
              if(is.null(dim(pairMat)) || ncol(pairMat)!=2) stop("pairMat must be matrix of 2 columns")
              if(!all(as.character(unique(c(pairMat[,1],pairMat[,2]))) %in% as.character(cl))) stop("Some elements of pairMat do not match cl")
              contrastNames <- apply(pairMat,1,function(y){y<-make.names(as.character(y));paste(y[1],y[2],sep="-")})
            }
            #     if(!removeNegative){
            #         levnames<-levels(cl)
            #         whNeg<-which(cluster<0)
            #         if(length(whNeg)>0){
            #             levnames[whNeg]<-paste("Neg",cluster[whNeg],sep="")
            #         }
            #     }
            #     if(removeNegative){
            #         levnames<-levels(factor(cluster[cluster>0]))
            #     }
            levnames<-make.names(levels(cl))
            if(outputType=="limma"){	
	            contr.matrix<-limma::makeContrasts(contrasts=contrastNames,levels=levnames)
				return(list(contrastMatrix=contr.matrix,contrastNames=names(contrastNames)))
			}
			if(outputType=="MAST") return(MAST::Hypothesis(contrastNames, levnames))
            
          })