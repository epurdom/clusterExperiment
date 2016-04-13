
#' Heatmap with two different sources for hierarchical clustering and color scale.
#' 
#' Make heatmap with color scale from one matrix and hiearchical clustering
#' from another. Also color palettes to go with heatmap
#' 
#' 
#' Note that plotHeatmap calles aheatmap under the hood. This allows you to
#' plot multiple heatmaps via par(mfrow=c(2,2)), etc. However, the dendrograms
#' do not resize if you change the size of your plot window in an interactive
#' session of R (this might be a problem for RStudio if you want to pop it out
#' into a large window...).
#' 
#' @name plotHeatmap
#' @aliases plotHeatmap
#' @docType methods
#' @details The plotHeatmap function calles \code{\link{aheatmap}} to draw the heatmap. The main point of \code{plotHeatmap} is to 1) allow for two different matrix inputs, one to visualize and one to cluster. 
#' 2) to assign colors to the clusters like in \code{\link{plotClusters}} that lines them up based on their similarity. 
#' The intended purpose is to allow the user to visualize the original count scale of the data (on the log-scale), but create the hierarchical clustering on another, more appropriate dataset for clustering, such as normalized data. Similarly, some of the palettes were developed assuming that the visualization might be on unscaled/uncentered data, rather than the residual from the mean of the gene, and thus palettes need to take on a greater range of relevant values so as to show meaningful comparisons with genes on very different scales.  
#'
#' @details If \code{annCol} contains a column of continuous data, whSampleDataCont should give the index of the column(s); otherwise the annotation data for those columns will be forced into a non-sensical factor (with nlevels equal the the number of samples). 
#'
#' @details If breaks is a numeric value between 0 and 1, then \code{breaks} is assumed to indicate the upper quantile (on the log scale) at which the heatmap color scale should stop. For example, if breaks=0.9, then the breaks will evenly spaced up until the 0.9 upper quantile of the log of the \code{heatData}, and then all values after the 0.9 quantile will be absorbed by the upper-most color bin. This can help to reduce the visual impact of a few highly expressed genes (variables). 
#'
#' @details If both \code{sampleData} and \code{clusterings} are provided, \code{clusterings} will be shown closest to data (i.e. on bottom) and \code{sampleData} on top. A user that wishes greater control should simply combine the two independently and give the combined matrix in sampleData.
#' @return Returns (invisibly) a list with elements that are passed to aheatmap.
#' \itemize{
#' \item{\code{breaks}}{The breaks used for aheatmap, after adjusting for quantile}
#' \item{\code{annCol}}{the annotation data.frame given to aheatmap}
#' \item{\code{clusterColors}}{the annotation colors given to aheatmap}
#' }
#' @author Elizabeth Purdom
#' @examples
#' 
#' data(simData)
#' data(simData)
#' cl<-rep(1:3,each=100)
#' cl2<-cl
#' changeAssign<-sample(1:length(cl),80)
#' cl2[changeAssign]<-sample(cl[changeAssign])
#' 
#' #simple, minimal, example. Show counts, but cluster on underlying means
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData)
#' 
#' #assign cluster colors
#' colors<-bigPalette[20:23]
#' names(colors)<-1:3
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,clusterColors=list(colors))
#' 
#' #show two different clusters
#' anno<-data.frame(cluster1=cl,cluster2=cl2)
#' out<-plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,annCol=anno)
#' #return the values to see format for giving colors to the annotations
#' out$clusterColors
#' 
#' #assign colors to the clusters based on plotClusters algorithm
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,annCol=anno,
#' alignClusterColors=TRUE)
#' 
#' #assign colors manually
#' annoColors<-list(cluster1=c("black","red","green"),
#' cluster2=c("blue","purple","yellow"))
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,annCol=anno,
#' clusterColors=annoColors)
#'
#' #give a continuous valued -- need to indicate columns
#' anno2<-cbind(anno,Cont=c(rnorm(100,0),rnorm(100,2),rnorm(100,3)))
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,annCol=anno2,
#' whSampleDataCont=3)

#' #compare changing breaks quantile on visual effect
#' \dontrun{
#' par(mfrow=c(2,2))
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,colorScale=seqPal1,
#' breaks=1,main="Full length")
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,colorScale=seqPal1,
#' breaks=.99,main="0.99 Quantile Upper Limit")
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,colorScale=seqPal1,
#' breaks=.95,main="0.95 Quantile Upper Limit")
#' plotHeatmap(cl,heatData=simCount,clusterSamplesData=simData,colorScale=seqPal1,
#' breaks=.90,main="0.90 Quantile Upper Limit")
#' }
#' 
#' 
#' 
setMethod(
    f = "plotHeatmap",
    signature = signature(data = "SummarizedExperiment"),
    definition = function(data, 
                          visualize=c("original","transformed","centeredAndScaled"),
                          orderSamples=c("hclust"), #can be indices too
                          whichFeatures=c("mostVar","all","PCA"), #can be indices too
                          nFeatures=NULL, sampleData=NULL, whSampleDataCont=NULL,
                          colorScale=if(centerAndScaleFeatures) seqPal3 else seqPal5,
                          isCount=FALSE,transData=NULL
    ){	
})
#colData of SE is data.frame, each row a sample, columns features.
#sampleData user gives either colnames(colData) or index of sample features; change sampleData to sampleData
#get rid of sampleData!
#whFeatures can be indices too, add gene names; also add lists of indices or lists of genes for having blanks between them.
setMethod(
  f = "plotHeatmap",
  signature = signature(data = "ClusterExperiment"),
  definition = function(data, 
                        visualize=c("original","transformed","centeredAndScaled"),
                       whichClusters= c("pipeline","all"),
                        orderSamples=c("dendrogramValue","hclust","orderSamplesValue","primaryCluster"), #can be indiecs too
                        whichFeatures=c("mostVar","all","PCA"), 
                       whSampleData=NULL,
                        nFeatures=NULL, 
                        colorScale=if(centerAndScaleFeatures) seqPal3 else seqPal5
  ){	
    visualize<-match.arg(visualize)
    
    .convertTry<-function(x,tryResult){if(!inherits(tryResult,"try-error")) return(tryResult) else return(x)}

    ####
    ##Set dimReduce, which will decide if 
    whichFeatures<-.convertTry(whichFeatures,try(match.arg(whichFeatures)))
    if(whichFeatures %in% c("mostVar","all","PCA")){ #
        dimReduce=switch(whichFeatures,
                         "mostVar"="mostVar",
                        "PCA"="PCA",
                        "all"="none")
        if(is.null(nFeatures)) nFeatures<-switch(whichFeatures,"mostVar"=500,"all"=NROW(sampleData),"PCA"=50)
        wh<-1:NROW(data)
    }
    else{
      if(is.character(whichFeatures)){#gene names
        if(is.null(rownames(data))) stop("Cannot give feature names in whichFeatures unless assay(data) has rownames")
        else{
          wh<-match(whichFeatures,rownames(data))
          if(all(is.na(wh))) stop("None of the feature names in whichFeatures match rownames(assay(data))")
          if(any(is.na(wh))){
            warning("Not all of the feature names in whichFeatures match rownames(assay(data))")
            wh<-na.omit(wh)
          }
        }  
      }
      else{#indices
        if(any(!whichFeatures %in% 1:NROW(data))) stop("invalid indices for whichFeatures")
        wh<-whichFeatures
      }
      dimReduce<-"none"
    }
    transObj<-.transData(transFun = transformation(assay(data[wh,]), nPCADims=nFeatures,nVarDims = nFeatures,dimReduce = dimReduce)
    if(dimReduce%in%"PCA") wh<-1:nFeatures
    if(dimReduce=="mostVar") wh<-transObj$whMostVar #give indices that will pull
    #########
    ##Assign heatData (i.e. visualizeData) 
    #########
    if(whichFeatures=="PCA") visualizeData<-transObj$x
    else{
      visualizeData<-switch(visualize,
                            "original"=assay(data[wh,]),
                            "transformed"=transObj$x,
                            "centeredAndScaled"=t(scale(t(transObj$x),center=TRUE,scale=TRUE))
                            )
    }
    
    ######
    #Create clusterSamplesData  
    ######
    orderSamples<-.convertTry(orderSamples,try(match.arg(orderSamples)))
    clusterSamples<-TRUE
    
    if(is.character(orderSamples)){
        if(orderSamples=="orderSamplesValues"){
            visualizeData<-visualizeData[,orderSamples(data)]
            clusterSamplesData<-visualizeData
            clusterSamples<-FALSE
        }
        if(orderSamples=="dendrogram"){
            if(length(dendrogram(sampleData))==0){
              ####Fix this once makeDendrogram available (#call clusterHclust on primary cluster)
                stop("No dendrogram assigned for this object")
                
            }
            else{
                clusterSamplesData<-dendrogram(sampleData)
            }
        }
        if(orderSamples=="hclust"){
            clusterSamplesData<-
        }
    }
    ##Determine which features:
    
    
  })
#assume that the clusterSamplesData,heatData has already been reduced in diminsionality in some way.

#' @param sampleData A matrix of additional sampleData to show with the heatmap. Unless indicated by \code{whSampleDataCont}, \code{sampleData} will be converted into factors.  ``-1'' indicates the sample was not assigned to a cluster and
#' gets color `unassignedColor' and '-2' gets the color 'missingColor'
#' @param whSampleDataCont Which of the sampleData columns are continuous and should not be converted to counts. NULL indicates no additional sampleData.
#' @param heatData matrix to define the color scale (e.g. the counts). Assumes samples on columns and variables on
#' rows.
#' @param clusterSamplesData Either a matrix that will be used to in hclust to define the hiearchical clustering of
#' samples (e.g. normalized data) or a pre-existing dendrogram for clustering the samples
#' @param clusterSamples Logical as to whether to do hierarchical clustering of
#' cells.
#' @param clusterSamples Logical as to whether to do any clustering of
#' features (otherwise will be in order given by heatData)
#' @param clusterFeatures Logical as to whether to do hiearchical clustering of
#' features.
#' @param showSampleNames Logical as to whether show cell names
#' @param showFeatureNames Logical as to whether show cell names
#' @param colorScale palette of colors for the color scale of heatmap
#' @param clusterColors Assignment of colors to the clusters. If NULL, clusters
#' will be assigned colors. If `annCol' should be list of length equal to
#' ncol(annCol) with names equal to the colnames of annCol; each element of the
#' list should be a vector of colors with names corresponding to the levels of
#' the column of annCol. If annCol=NULL (i.e. use clusterVec), then the name of
#' the length-1 list should be `Cluster'
#' @param alignClusterColors Logical as to whether should align the clusters when
#' colors are assigned (only used if annCol=NULL)
#' @param breaks Either a vector of breaks (should be equal to length 52), or a
#' number between 0 and 1, indicating that the breaks should be equally spaced
#' (based on the range in the heatData) upto the `breaks' quantile, see details.
#' @param unassignedColor color assigned to cluster values of '-1' ("unassigned")
#' @param missingColor color assigned to cluster values of '-2' ("missing")
#' @param ... passed to aheatmap
#' @rdname plotHeatmap
setMethod(
    f = "plotHeatmap",
    signature = signature(data = "matrix"),
    definition = function(data,sampleData, 
                          clusterSamplesData=heatData, 
                          clusterFeaturesData=heatData, 
                          whSampleDataCont=NULL,
                          clusterSamples=TRUE,showSampleNames=FALSE, 
                          clusterFeatures=TRUE,showFeatureNames=FALSE,
                          colorScale=seqPal5,
                          clusterColors=NULL,alignClusterColors=FALSE,
                          unassignedColor="white",missingColor="grey", breaks=NA,...
    ){	
        
        ##########
        ##Deal with clusterings...
        ##########
        #check clusterings input:
        if(!is.matrix(sampleData) | !is.data.frame(sampleData)){ stop("sampleData must be a either a matrix or a data.frame") 
        }
        if(NCOL(heatData) != NROW(sampleData)) stop("sampleData must have same number of rows as columns of heatData")	
        if(!is.null(sampleData)){
            ###Make sampleData explicitly factors, except for whSampleDataCont 
            ###(not sure why this simpler code doesn't give back data.frame with factors: annCol<-apply(annCol,2,function(x){factor(x)}))
            
            tmpDf<-do.call("data.frame",lapply(1:ncol(sampleData),function(ii){factor(sampleData[,ii])}))
            names(tmpDf)<-names(sampleData)
            if(!is.null(whSampleDataCont)) tmpDf[,whSampleDataCont]<-sampleData[,whSampleDataCont]
            annCol<-tmpDf
            
            if(is.null(clusterColors)){ #assign default colors
                if(is.null(whAnnCont) || length(whAnnCont)<ncol(annCol)){ #if not all continuous
                    #Pull out the factors to assign clusters
                    if(!is.null(whSampleDataCont)) tmpDf<- annCol[,-whAnnCont,drop=FALSE] else tmpDf<-annCol
                    if(alignClusterColors){
                        #align the clusters and give them colors
                        alignObj<-plotClusters(data.matrix(tmpDf) ,plot=FALSE,unassignedColor=unassignedColor, missingColor=missingColor) 
                        annColors<-convertClusterColors(alignObj$clusterColors,output=c("aheatmapFormat"))
                    }
                    else{#give each distinct colors, compared to row before
                        maxPerAnn<-sapply(tmpDf,function(x){max(as.numeric(x))})
                        annColors<-mapply(tmpDf,c(0,head(cumsum(maxPerAnn),-1)),FUN=function(fac,add){
                            cols<-.thisPal[1:nlevels(fac)+add]
                            names(cols)<-levels(fac)
                            cols[names(cols)=="-1"]<-unassignedColor #unassigned 
                            cols[names(cols)=="-2"]<-missingColor #
                            cols<-cols[order(names(cols))]
                            return(cols)
                        },SIMPLIFY=FALSE)
                    }
                    
                }
            }
            else{
                #check if not in our format, assume in aheatmap format
                if(!any(sapply(clusterColors,function(x){is.null(dim(x))}))) annColors<-convertClusterColors(clusterColors,output=c("aheatmapFormat"))
                else annColors<-clusterColors
            }
        }
        else annCol<-NULL
        
        ##########
        ##Deal with numeric matrix for heatmap ...
        ##########
        heatData<-data.matrix(heatData) 
        ###Create the clustering dendrogram:
        if(clusterSamples){
            if(inherits(clusterSamplesData, "dendrogram")){
                if("nobs.dendrogram" %in% methods(nobs)){ #not all versions have nobs method for dendrograms; R3.2.0 has it; 3.1.1 doesn't
                    if(nobs(clusterSamplesData)!=ncol(heatData)) stop("clusterSamplesData dendrogram is not on same number of observations as heatData")
                }
                else{warning("This R version doesn't doesn't allow for checking that dendrogram supplied has correct number observations. If not, surprising downstream errors can occur.")}
                dendroCells<-clusterSamplesData	
            } 
            else{
                clusterSamplesData<-data.matrix(clusterSamplesData)
                #check valid
                if(ncol(clusterSamplesData)!=ncol(heatData)) stop("clusterSamplesData matrix does not have on same number of observations as heatData")
                dendroSamples<-as.dendrogram(hclust(dist(t(clusterSamplesData)))) #dist finds distances between rows
            }
        }
        else{ 
            clusterSamples<-NA
        }
        if(clusterFeatures){
            if(inherits(clusterFeaturesData, "dendrogram")){
                if("nobs.dendrogram" %in% methods(nobs)){ #not all versions have nobs method for dendrograms; R3.2.0 has it; 3.1.1 doesn't
                    if(nobs(clusterFeaturesData)!=ncol(heatData)) stop("clusterSamplesData dendrogram is not on same number of observations as heatData")
                }
                else{warning("This R version doesn't doesn't allow for checking that dendrogram supplied has correct number observations. If not, surprising downstream errors can occur.")}
                dendroFeatures<-clusterFeaturesData	
            } 
            else{
                clusterFeaturesData<-data.matrix(clusterFeaturesData)
                #check valid
                if(ncol(clusterFeaturesData)!=ncol(heatData)) stop("clusterFeaturesData matrix not have on same number of observations as heatData")
                dendroFeatures<-as.dendrogram(hclust(dist(clusterFeaturesData))) #dist finds distances between rows
            }
        }
        else{ 
            clusterFeatures<-NA
        }
        #############
        # put into aheatmap
        #############
        breaks<-setBreaks(breaks,heatData)
        out<-NMF::aheatmap(heatData,  
                           Rowv =if(!is.na(clusterFeatures) && clusterFeatures) dendroFeatures else clusterFeatures,  
                           Colv = if(!is.na(clusterSamples) && clusterSamples) dendroSamples else clusterSamples, 
                           color = colorScale, scale = "none",
                           annCol = annCol,annColors=annColors,breaks=breaks,...)
        
        #############
        # add labels to clusters at top of heatmap
        #############
        #	browser()
        if(!is.null(annCol)){
            newName<-NMF:::vplayout(NULL) #will be 1 greater (hopefully!) this is fragile. Don't know if it will always work.
            newNameList<-strsplit(newName,"\\.")[[1]]
            oldIndex<-as.numeric(newNameList[[3]])-1
            newNameList[[3]]<-oldIndex
            oldName<-paste(newNameList,collapse=".")
            grid::seekViewport(sprintf("aheatmap-%s",oldName))
            NMF:::vplayout(3,4:5)
            #grid::grid.rect()
            y <- seq(0,1,length=ncol(annCol))
            n<-ncol(annCol)
            y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
            #		grid::grid.points(x = grid::unit(rep(0,length(y)),"npc"),y = grid::unit(y[n:1], "bigpts"))
            grid::grid.text(colnames(annCol), x = grid::unit(rep(0.05,length(y)),"npc"),y = grid::unit(y[n:1], "bigpts"), vjust = 0.5, hjust = 0,gp= grid::gpar(fontsize=10))
            grid::upViewport() #close it
            grid::upViewport() #close it
        }
        
        invisible(list(aheatmapOut=out,sampleData=annCol,clusterColors=clusterColors,breaks=breaks))
    }
)
