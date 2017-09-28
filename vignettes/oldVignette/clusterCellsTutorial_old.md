<!-- ---
title: "clusterExperiment Tutorial"
author: "Calvin Chi and Elizabeth Purdom"
date: "2016-03-04"
output: html_document
--- -->

# Introduction
The goal of this package is to allow the user to try many different clustering algorithms in one package structure, while implementing common post-processing steps unrelated to the clustering algorithm (e.g. consensus clustering). The package also provides tools for doing differential expression analysis for the clusters found to find important genes, as well as visualization techiniques. 



## Basic overview of clustering routines

The package encodes many common practices that are shared across clustering algorithms, like subsampling the data, not clustering samples with negative silhouete scores, sequentially removing clusters and reclustering, and so forth, and allows the user to simply make different choices the parameters as well as the underlying clustering algorithm. 

There are two main user-functions for clustering, `clusterAll` and `compareChoices`. `clusterAll` is the wrapper function that calls the underlying clustering functions. `compareChoices` is a convenience function that implements `clusterAll` across combinations of parameter choices. Such parameter choices can be the clustering method, whether to subsample the data, or the specific parameters of the clustering algorithm (e.g. the number of clusters in kmeans).

A common workflow could be to apply `compareChoices` to the data to get a large collection of clusterings, each from clustering based on different parameters. To try to find a clustering presumably robust to the parameter choices, `clusterExperiment` provides a function `findSharedClusters` that tries to find a unifying clustering of the samples. 

We find that many methods for choosing the appropriate number of clusters for methods like k-means err on the side of smaller number of clusters. However, we find in practice that we tend to prefer to err on finding many clusters and then merging them based on examining the data (and many of the default parameters are set for finding many clusters). We provide the function `clusterHclust` to perform hierarchical clustering of the resulting clusters and the function `mergeClusters` to follow that hiearchy and merge the clusters based on the percentage of significantly different genes found in the comparisons of different arms of the hiearchy.


## Finding related genes

A common practice after determining clusters is to perform differential gene expression analysis between the clusters in order to find genes that show the greatest differences amongst the clusters. We would stress that this is purely an exploratory technique, and any p-values that result from this analysis are invalid, in the sense that they are likely to be inflated. This is because the same data was used to define the clusters as to perform differential expression analysis. 

Since this is a common task, we provide the function `getBestGenes` to perform various kinds of differential expression analysis between the clusters. A common F-statistic is performed. However, we find that it is far more informative to do pairwise comparisons between clusters, or one cluster against all in order to find genes that are specific to a particular cluster. An option for all of these choices is provided in the `getBestGenes` function. 

In addition, the `getBestGenes` function provides the ability to do a "voom" correction to account for the mean-variance relationship that is common in count data (or FPKM or TPM). Unlike edgeR or DESeq, the voom correction does not require a count matrix, and therefore can be used on FPKM or TPM entries, as well as normalized data. 

## Visualization


## Data example

We will work with a simulated gene expression data set with 300 samples and 3 hypothetical genes, which we will simply label as *growth*, *immune*, and *structural* genes. The expression values will be drawn from different normal distributions. Hence, our data matrix has dimensions $300 x 3$. In this data set, we know *a priori* that there are k = 5 clusters (1-100, 101-125, 126-175, 176-225, 226-300). There are 100 samples in the first cluster, 25 in the second cluster, 50 samples in the third cluster, 50 samples in the fourth cluster, and 75 samples in the fifth cluster. Let us first visualize the clusters. 


```
## Error in library(scatterplot3d): there is no package called 'scatterplot3d'
```

```
## Error in eval(expr, envir, enclos): could not find function "scatterplot3d"
```

```
## Error in strwidth(legend, units = "user", cex = cex, font = text.font): plot.new has not been called yet
```


# `ClusterAll`

Let us start with a walkthrough of `clusterAll`. As a wrapper function, `clusterAll` can receive a variety of arguments to perform different clustering algorithms. (The options are documented in the more detailed documentation vignette). Here are the argument calls for `clusterAll`: 

```r
clusterAll(x, subsample = TRUE, sequential = FALSE, clusterFunction = c("tight", "hierarchical01", "pam", "kmeans"), clusterDArgs = NULL, subsampleArgs = NULL, seqArgs = NULL)
```
These arguments can be better organized and visualized as a diagram: 
<img src="R_figure/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

## PAM Clustering
Let us use the function `clusterAll` from the `clusterExperiment` package to perform simple `pam` clustering on our expression data set with `k=5`. In this expression matrix, rows are samples and columns are features. In the first example of pam clustering, `subsample=FALSE` and `sequential=FALSE`. If `subsample=FALSE`, then `clusterFunction` must be `pam`. pam clustering with `clusterAll` returns a one-element list with a vector of cluster assignments for each sample.
<img src="R_figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```r
library(clusterExperiment)
library(cluster)
load("data/expressions.Rda")
```

```
## Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```

```r
expression = expressions[, 1:3]
simpleCluster = clusterAll(expression, subsample=FALSE, sequential=FALSE, clusterFunction="pam", clusterDArgs=list('k'=5))
table(simpleCluster$clustering)
```

```
## 
##   1   2   3   4   5 
## 100  75  50  50  25
```
Let us compare this result with the result from calling `pam` with `k=5` from the `cluster` package. 

```r
pamCluster = pam(dist(expression), 5)
tab = table(simpleCluster$clustering, pamCluster$clustering)
barplot(tab, xlab="Clusters", ylab="Counts", main="clusterAll vs pam", legend=rownames(tab))
```

<img src="R_figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" style="display: block; margin: auto;" />
We have now verified that these two function calls perform the same task. 

Now let us perform pam clustering with adjustments to the parameters `findBestK`, `kRange`, and `removeSil`. These parameters are passed to the argument `clusterDArgs` in `clusterAll`. The default cutoff for the minimum silhouette width is 0. This cutoff can be adjusted via the parameter `silCutoff` passed to `clusterDArgs`. 

```r
Cluster<-clusterAll(expression, subsample=FALSE, sequential=FALSE, clusterFunction="pam", clusterDArgs=list(findBestK=TRUE, removeSil=TRUE, kRange=2:10))

table(Cluster$clustering)
```

```
## 
##  -1   1   2   3 
##   3 154  93  50
```
## Subsampling
When `subsample=TRUE`, `clusterAll` will first subsample the rows of the data matrix, cluster the subsamples using a certain clustering algorithm, and return a clustering co-occurance matrix. Then, `clusterAll` will cluster the co-occurance matrix using another specified algorithm, which can be the same as the algorithm used in subsampling. 

In this demonstration let us cluster the subsamples using `kmeans` with `k=5`, then cluster the co-occurance matrix using the tight clustering algorithm. 

<img src="R_figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" />

```r
Cluster<-clusterAll(expression, subsample=TRUE, sequential=FALSE, clusterFunction="tight", subsampleArgs=list("k"=5, clusterFunction="kmeans"))

table(Cluster$clustering)
```

```
## 
## -1  1  2  3  4  5  6  7  8 
## 18 65 60 47 39 25 25 18  3
```

## Sequential 
When `sequential=TRUE`, `clusterExperiment` performs sequential clustering, namely clustering the data, removing the cluster, and reclustering the remaining data. 

In this demonstration of sequential clustering, the necessary argument `k0=5` is passed to `seqArgs` as the value of K at the first iteration of sequential algorithm. 
<img src="R_figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" style="display: block; margin: auto;" />

```r
Cluster<-clusterAll(expression, subsample=TRUE, sequential=TRUE, clusterFunction="hierarchical01", seqArgs=list(k0=5, verbose = FALSE), subsampleArgs = list(clusterFunction="kmeans"))
```
Running `clusterAll` with `sequential=TRUE` will return a list of the three elements `clustering`, `clusterInfo`, and `whyStop`. `clusterInfo` is a matrix of information regarding the algorithm behavior for each cluster (the starting and stopping K for each cluster, and the number of iterations for each cluster). `whyStop` is a character string explaining what triggered the algorithm to stop.

```r
table(Cluster$clustering)
```

```
## 
##  -1   1   2   3   4   5   6   7   8   9  10  11  12  13 
##  28  23  27   1   1  67 100   8   1   5  11  11  14   3
```

```r
Cluster$clusteringInfo
```

```
## NULL
```

```r
Cluster$whyStop
```

```
## [1] "Ran out of samples"
```
******
# Processing clusterAll Output with getBestGenes
The function `getBestGenes` calls limma on input data to determine the gene features most associated with found clusters. The minimal arguments to `getBestGenes` are `type` and `contrastAdj`, which indicate what tests to perform and what type of FDR correction to do for contrast tests respectively. 

In this demonstration, `getBestGenes` will be called to determine the features most associated with the clusters found in a new 300x50 dataset called `simData`. Clustering will first be performed with PAM with `k=4`. 

```r
data(simData)
cl = clusterAll(simData, clusterFunction="pam", subsample=FALSE, sequential=FALSE, clusterDArgs=list(k=4))
pairsAll<-getBestGenes(cl$clustering,simData,type="Pairs",contrastAdj="All")

head(pairsAll)
```

```
##   IndexInOriginal Contrast  Gene      logFC     AveExpr         t
## 1              20    X1-X2 Row20 -10.756333  0.02714648 -22.65527
## 2              30    X1-X2 Row30 -10.105420  0.32693884 -22.59699
## 3               5    X1-X2  Row5 -10.882924  0.23130864 -21.77647
## 4              32    X1-X2 Row32 -10.237611 -0.00352073 -21.71944
## 5               2    X1-X2  Row2 -10.114661  0.20750472 -21.08549
## 6              29    X1-X2 Row29  -9.779938  0.05160590 -20.90595
##        P.Value    adj.P.Val        B
## 1 8.621062e-67 6.435669e-64 141.5543
## 2 1.402107e-66 6.435669e-64 141.0698
## 3 1.366599e-63 4.181792e-61 134.2127
## 4 2.210145e-63 5.072283e-61 133.7337
## 5 4.710809e-61 8.649045e-59 128.3908
## 6 2.164024e-60 3.310957e-58 126.8715
```
The significantly associated features are indicated in the `IndexInOriginal` column of the dataframe returned by `getBestGenes`. One way to visualize samples with only the features in `IndexInOriginal` is to use the `dualHeatmap()` function. The `dualHeatmap()` function minimally requires two arguments - a vector of cluster assignments and the data matrix. `dualHeatmap()` is different from regular heatmaps in that it also plots cluster assignments on the top, where different colors denote different cluster assignments. In the code example below, `clusterData` is set to the data matrix to define the color scale. 

```r
dataMatrix = simData[, unique(pairsAll$IndexInOriginal)]
cluster = cl$clustering
dualHeatmap(cluster, dataMatrix, clusterData=dataMatrix)
```

<img src="R_figure/unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" style="display: block; margin: auto;" />

```
## 
## The downloaded binary packages are in
## 	/var/folders/h4/xtpbyfq55qd3rc882bm4zfjw0000gn/T//RtmpDKths0/downloaded_packages
```


******
# compareChoices Workflow

`compareChoices` is a function for running many different combinations of parameters or different datasets. Before running clustering on these combinations, it may be a good idea to inspect all the parameter combinations by setting `run`=FALSE. In this example, the combination of parameters we are testing are the range of k values for clustering and whether `removeSil` is TRUE or FALSE. 

```r
parameters = compareChoices(expression, ks=c(3:7), clusterMethod=c("pam", "kmeans"), run=FALSE)
```

```
## 10 parameter combinations, 0 use sequential method.
```

```r
head(parameters)
```

```
##                                     dataset k alpha findBestK sequential
## k=3,alpha=NA,clusterMethod=pam     dataset1 3    NA     FALSE      FALSE
## k=4,alpha=NA,clusterMethod=pam     dataset1 4    NA     FALSE      FALSE
## k=5,alpha=NA,clusterMethod=pam     dataset1 5    NA     FALSE      FALSE
## k=6,alpha=NA,clusterMethod=pam     dataset1 6    NA     FALSE      FALSE
## k=7,alpha=NA,clusterMethod=pam     dataset1 7    NA     FALSE      FALSE
## k=3,alpha=0.1,clusterMethod=kmeans dataset1 3   0.1     FALSE      FALSE
##                                    removeSil subsample clusterMethod
## k=3,alpha=NA,clusterMethod=pam         FALSE     FALSE           pam
## k=4,alpha=NA,clusterMethod=pam         FALSE     FALSE           pam
## k=5,alpha=NA,clusterMethod=pam         FALSE     FALSE           pam
## k=6,alpha=NA,clusterMethod=pam         FALSE     FALSE           pam
## k=7,alpha=NA,clusterMethod=pam         FALSE     FALSE           pam
## k=3,alpha=0.1,clusterMethod=kmeans     FALSE     FALSE        kmeans
##                                    silCutoff
## k=3,alpha=NA,clusterMethod=pam             0
## k=4,alpha=NA,clusterMethod=pam             0
## k=5,alpha=NA,clusterMethod=pam             0
## k=6,alpha=NA,clusterMethod=pam             0
## k=7,alpha=NA,clusterMethod=pam             0
## k=3,alpha=0.1,clusterMethod=kmeans         0
```
Let us run `compareChoices` on the original `expression` data set we simulated.

```r
result = compareChoices(expression, ks=c(3:7), clusterMethod=c("pam"), subsample=FALSE, removeSil = c(TRUE, FALSE))
```

```
## 10 parameter combinations, 0 use sequential method.
```

```r
head(result$clMat)
```

```
##      k=3,removeSil=1 k=4,removeSil=1 k=5,removeSil=1 k=6,removeSil=1
## [1,]               1               1               1               2
## [2,]               1               1               1               2
## [3,]               1               1               1               2
## [4,]               1               1               1               5
## [5,]               1               1               1               5
## [6,]               1               1               1               2
##      k=7,removeSil=1 k=3,removeSil=0 k=4,removeSil=0 k=5,removeSil=0
## [1,]               3               1               1               1
## [2,]               3               1               1               1
## [3,]               3               1               1               1
## [4,]               1               1               1               1
## [5,]               1               1               1               1
## [6,]               3               1               1               1
##      k=6,removeSil=0 k=7,removeSil=0
## [1,]               2               2
## [2,]               2               2
## [3,]               2               2
## [4,]               4               1
## [5,]               4               1
## [6,]               2               2
```
Actually running `compareChoices` will return a list with elements `clMat` and `clusterInfo`. `clusterInfo` will be a list with information regarding clustering result only when sequential=TRUE while `clMat` will be a matrix with each row representing a sample and each column the cluster assignments from a particular combination of parameters. 

After `compareChoices` is run, the function `plotTracking` can be used to visualize a color-coded plot of cluster assignments across different combinations of parameters. `plotTracking` minimally takes as its argument a matrix of cluster assignments where each row represents a sample and each column a clustering. 
<img src="R_figure/unnamed-chunk-18-1.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" style="display: block; margin: auto;" />

```r
par(mar=c(2, 8, 1, 1))
plotTracking(result$clMat)
```

<img src="R_figure/unnamed-chunk-19-1.png" title="plot of chunk unnamed-chunk-19" alt="plot of chunk unnamed-chunk-19" style="display: block; margin: auto;" />

```r
par(mar=c(4, 4, 1, 1))
```

To find a consensus cluster across many different clusterings of the same data `findSharedClusters` can be used. The minimal input to `findSharedClusters` is a cluster matrix with samples as rows and different clustering assignments for each column. 
<img src="R_figure/unnamed-chunk-20-1.png" title="plot of chunk unnamed-chunk-20" alt="plot of chunk unnamed-chunk-20" style="display: block; margin: auto;" />
As an example, let us find the consensus cluster from the output of `compareChoices`. 

```r
r = findSharedClusters(result$clMat)
table(r)
```

```
## Error in base::table(...): all arguments must have the same length
```
******
#Change Arguments
The functions in package `clusterAll` accept a myraid of arguments, and while most parameters have default settings, it is important to know how to change arguments to suit a particular need. The best way to learn about additional arguments for a function is to read the appropriate help files by executing `?` followed by the name of the function. 

Let us look at an example of changing arguments in `clusterAll` with `clusterFunction="tight"`, `subsample=TRUE`, and `sequential=FALSE`. 

```r
cluster<-clusterAll(expression, subsample=TRUE, sequential=FALSE, clusterFunction="tight", clusterDArgs=list(alpha=0.05), subsampleArgs=list("k"=5, clusterFunction="kmeans"))
table(cluster$clustering)
```

```
## 
## -1  1  2  3  4  5  6  7  8  9 10 11 
## 20 49 47 46 43 28 21 18 17  6  3  2
```
In the code just executed, the arguments fed to `clusterAll` are minimum arguments to make this particular clustering algorithm run. However, the rest of the arguments are set at default value. 

Pretend that the user is interested in adjusting the subsampling process so that the proportion of subsamples is 0.3 and the number of times sampled is 10. How would he or she learn how to change the arguments in `clusteAll` to execute the desired change? A good way to begin would be to access the help file for `clusterAll`.

```r
?clusterAll
```
A help file will have many sections. The `Usage` section will illustrate how to pass arguments into a function and indicate any default values. The `Arguments` section will explain what each argument means. In our example, because we are adjusting the subsampling procedure, feeding arguments to `subsampleArgs` seems like a reasonable choice. The `Arguments` section of the help file explains that `subsampleArgs` is: 

> list of arguments to be passed to subsampleClustering

What are those arguments? It appears that `clusterAll` will pass arguments in `subsampleArgs` to the function `subsampleClustering` during execution, so it will be informative to look at the help file for `subsampleClustering`. 

```r
?subsampleClustering
```
Now the help file for `subsampleClustering` explains that the arguments `resamp.num` and `samp.p` are the arguments that set the number of subsamples to draw and the proportion of samples to draw respectively. We now know what arguments to pass to the original `subsampleArgs` in `clusterAll`. 

```r
cluster<-clusterAll(expression, subsample=TRUE, sequential=FALSE, clusterFunction="tight", clusterDArgs=list(alpha=0.05), subsampleArgs=list("k"=5, clusterFunction="kmeans", resamp.num=10, samp.p=0.3))
table(cluster$clustering)
```

```
## 
## -1  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
## 11 56 46 36 30 30 25 20 12  9  7  5  3  3  3  2  2
```

Learning how to changing other arguments for other functions will likely invoke a similar logic, and it is up to the user to decide what to explore. 


******
# Generate co-occurance matrix
To generate the co-occurance matrix after clustering subsamples of the data matrix, one can call the `subsampleClustering` function. The main arguments for the `subsampleClustering` are as follows: 

`subsampleClustering` returns a $n x n$ matrix of probability of co-occurance in the same cluster, with $n$ = number of samples. Let us generate a co-occurance matrix from our gene expression data set with 100 subsamples, each composed of 70% of the data. 


```r
require(gplots)
```

```
## Loading required package: gplots
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
d = subsampleClustering(expression, k=5, clusterFunction="kmeans", resamp.num=100, sam.p=0.7)
heatmap.2(d, col=redgreen(75), scale="none", density.info="none", trace="none") 
```

<img src="R_figure/unnamed-chunk-26-1.png" title="plot of chunk unnamed-chunk-26" alt="plot of chunk unnamed-chunk-26" style="display: block; margin: auto;" />
