0.0.0.9006:

* fixed so that mergeClusters, clusterHclust, and getBestGenes will appropriately convert if the input of clustering vector is a factor rather than numeric (with warning).
* fixed mergeClusters to have option to indicate that input matrix is a count matrix (in which case will create dendrogram with log(counts+1) and will do getBestGenes with the voom correction) 
* added more tutoral-oriented vignette (old vignette is now the documentation vignette with more detail about the internal workings of package). Currently is just simulated data, but will be updated to real single-cell sequencing dataset.

0.0.0.9005

* Changed simulated data so load all with data(simData) rather than separate calls for simData and simCount. Also added 'trueCluster' vector to give true cluster assignments of simulated data
* added dendro example to getBestGenes
* added example to clusterHclust
* added single function for converting to phylobase tree (used internally by package)
* added functionality to find proportion of significant null hypotheses for merging clusters (mergeClusters)

0.0.0.9004:

* Changed compareChoices.R to only set k<-NA if sequential=FALSE (previously for all where findBestK=TRUE)
* Added to vignette
* fixed bug in plotTracking to correctly plot "-1"