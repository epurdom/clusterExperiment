Unreleased:
* Changed simulated data so load all with data(simData) rather than separate calls for simData and simCount. Also added 'trueCluster' vector to give true cluster assignments of simulated data
* added dendro example to getBestGenes
* added single function for converting to phylobase tree
* added functionality to find proportion of significant null hypotheses for merging clusters

0.0.0.9004:

* Changed compareChoices.R to only set k<-NA if sequential=FALSE (previously for all where findBestK=TRUE)
* Added to vignette
* fixed bug in plotTracking to correctly plot "-1"