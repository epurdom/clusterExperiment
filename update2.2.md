# New to Version 2.2

**Important** If you have objects created with 2.0.0, you should run `updateObject` to update the class definition because there have been changes to the class definition since 2.0.0:

```
ceObj<-updateObject(ceObj)
```


There have also been a number of changes and enhancements to the package. These are the most important

* We have changed the function `combineMany` to `makeConsensus`. This has resulted in changes to the names of the arguments of `RSEC`
	- `combineProportion` -> `consensusProportion` in `RSEC`
	- `combineMinSize` -> `consensusMinSize` in `RSEC`
* Add functionality to `getBestFeatures` to allow `edgeR` for DE, as well as weights used with `edgeR` for compatability with weights to handle zero-inflation. As part of this change  `isCount` argument has been replaced with more fine-grained `DEMethod` argument in `getBestFeatures`, `mergeClusters`; and the argument `mergeDEMethod` in `RSEC` is now available.
* We have changed the argument `sampleData` in various plotting commands to `colData` to better indicate that the argument is to identify columns in `colData` that should also be plotted. Furthermore `plotDendrogram` now takes the argument `colData` for plotting of information in `colData` with the dendrogram.
* We have changed the names of arguments related to unassigned (`-1` or `-2` assignments) to more consistently use the term "unassigned", as well as adding the function `assignUnassigned`:
	- argument `removeNegative` -> `removeUnassigned` in `getBestFeatures` 
	- argument `ignoreUnassignedVar` -> `filterIgnoresUnassigned` in `mergeClusters` (and other functions) for clarity.
	- function `removeUnclustered` -> `removeUnassigned`
* New plotting functions:
	- `plotTableClusters`
	- `plotFeatureScatter`
* Allow the arguments `subsample` and `sequential` to `RSEC` to allow for opting out of those options for large datasets (but default is `TRUE` unlike `clusterMany`)
* The argument `whichAssay` is added to most functions to allow the user to select the assay on which the operations will be performed.


These changes are fully detailed in the [NEWS](https://github.com/epurdom/clusterExperiment/blob/master/NEWS) file of the package (all releases since May 1). 




