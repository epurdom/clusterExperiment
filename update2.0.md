# New to Version 2.0

In version 2.0.0 of `clusterExperiment` we have made some major changes to the package. The most relevant changes for the user are:

* Changed how we handle dimensionality reduction and filtering of genes -- this includes changes to the names of important arguments of key functions!
* Added support for the main data matrix to be stored in a HDF5 file (i.e. not read into memory)
* Improved implementation of subsampling (in C++) to improve speed 
* Some functions changed names 
* Some additional helper functions were added

These changes are fully detailed in the [NEWS](https://github.com/epurdom/clusterExperiment/blob/master/NEWS) file of the package (all releases since October 31, 2017). 

Below we will explain the first two of these changes in more detail. We would also note that in Bioconductor 3.6 (the release before) there were many more plotting functions that were introduced. 

## Dimensionality reduction and filtering

### Overview

This version consists of a major update of how dimensionality reduction and filtering is done. The `ClusterExperiment` class has been updated to extend the new bioconductor class, `SingleCellExperiment` (see the [SingleCellExperiment package](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)). 

Briefly, the `SingleCellExperiment` class extends the `SummarizedExperiment` class to give a structure for saving the results of dimensionality reductions.  This is done by adding a slot `reducedDims`, which is a `SimpleList` of datasets that have the same number of observations but reduced dimensions (i.e. genes). This gives a unified slot for saving the results of applying a dimensionality reduction method and helper functions to access them, etc. The user gives them names, e.g. "PCA" or "tSNE". 

**Dimensionality Reduction** `clusterExperiment` now makes use these stored dimensionality reductions in functions like `RSEC` and `clusterMany`. This allows `clusterExperiment` to make use of any dimensionality reduction method so long as the user saves it in the appropriate slot in a `SingleCellExperiment` object. The user can also choose not to calculate any dimensionality reduction and just do as before and have a function like `clusterMany` do the dimensionality reduction (i.e. PCA) internally. The difference is that now the results of the PCA will be stored in the appropriate slot so that they will not need to be recalculated in the future. 

**Filtering of genes** We also added in `clusterExperiment` package a similar option for filtering statistics. In particular, in `clusterExperiment` we have always allowed the user -- instead of dimensionality reduction via PCA-like methods -- to instead reduce the dimensionality of the problem by filtering to the top set of genes, e.g. the top 500 most variable genes. In this case `clusterExperiment` will calculate the variance for every gene and reduce down to the top 500 genes. Now `clusterExperiment` when calculating statistics (like `var` or `mad`) will add the per-gene value of the statistic in a column of the `rowData` of the `ClusterExperiment` object (`rowData` is a standard slot of a `SummarizedExperiment`). Similarly, if the user has already calculated a per-gene statistic and saved it as a column in the `rowData` slot, this user-defined statistic can be used for filtering. This means that the user is not limited to the built-in functions provided in `clusterExperiment`.

Note that simplicity we call "dimensionality reduction methods" to be those like PCA that reduce the data in a way that is *not* a simple selection of existing variables, but rather create new variables to represent the data. This is because a simple selection of variables can be stored as a single vector of the length of the number of genes and the reduced data can be obtained from the original matrix. The more complicated methods actually have to save a matrix with a value for each observation for each new variable. 

### Details 

A great deal has changed under the hood of any function that allowed for dimensionality reduction and greatly simplified and unified our treatment of filtering and dimensionality reduction. The main function affected is `clusterMany`, which runs the clustering, but  `makeDendrogram` is another example (and of course `RSEC` which is a wrapper around these).

To make the function compatible with `SingleCellExperiment`, we have changed many of the names of our arguments related to dimensionality reduction. This is because the slot names and related functions of `SingleCellExperiment` take the form of `reducedDims` and our previous versions of `clusterExperiment` used `dimReduce` format instead. We also changed to the name of arguments to be less tied to 'PCA' and 'var':

- `nPCADims` changed to `nReducedDims` in clusterMany-related functions
- `nVarDims` changed to `nFilterDims` in clusterMany-related functions
- `dimReduce` changed to `reduceMethod` across functions
- `plotDimReduce` to `plotReducedDims`
- `ndims` to `nDims` in `clusterSingle` and `makeDendrogram` to keep consistency.

The package `clusterExperiment` has built in functions for both dimensionality reduction (right now only PCA) and for filtering. These can be obtained by the new functions `listBuiltInReducedDims` and `listBuiltInFilterStats`, which give a character vector of the names of currently available functions for dimensionality reduction and statistics for filtering genes, respectively.

The argument `reduceMethod`, like the previous argument `dimReduce`, defines *either* the dimensionality reduction method or the filtering method to be used to reduce the number of dimensions of the data that will be used (and for `clusterMany` this can be a mixture of the two if the user wants to compare them). In places where multiple values can be given  (i.e. `clusterMany` or `RSEC`) the user has the choice to give to `reduceMethod` either 

* the names of  stored values (which can be a mixture of names of stored statistics for filtering and names of stored dimensionality reduction methods) OR 
* names of built-in functions provided by `clusterExperiment` to be calculated internally (as given by `listBuiltInReducedDims` and `listBuiltInFilterStats`)

The user cannot do both (i.e. give `reducedMethod` 2 names that match user-defined stored values and 3 names that are built-in functions). To do this the user can call the new functions `makeReducedDims` and `makeFilterStats` that will apply the built-in method (and store them appropriately) the built-in methods for dimensionality reductions and statistics for filtering the data, respectively. These are indeed the functions called internally by functions like `clusterMany`. In the above example, after calling these functions for the 3 built-in functions, the user can then call `clusterMany` on all 5 of the names that are stored values (the 2 user defined, and the 3 that were created by the built-in functions).


## HDF5 Support

The `clusterExperiment` packages is now compatible with `SummarizedExperiment` objects that have `DelayedArray` classes in their `assay` slot (which includes `HDF5Matrix` and  `DelayedMatrix`). These are classes from the `DelayedArray` and `HDF5Array` packages that allow the assay to be stored on file rather than in memory. 

Note, however, that while the package allows for these objects, it doesn't mean that it makes use of the HDF5 structure. Many times if the `assay` object must be actually used beyond simple subsetting, it will call the entire matrix into memory for the computations. The advantages currently of having the full dataset in HDF5 format are:

1. If you are not using the full assay for calculations, but rather a version after dimensionality reduction, then the full assay will not be needed in memory, but only the reduced version of the data. This is true for filtering as well as the more complicated dimensionality reductions stored in the `reducedDims` slot, since subsetting of the HDF5 matrix doesn't call the entire matrix into memory. Furthermore, *some* (but not all) of our built-in filtering statistics are HDF5 aware, and do not call the entire matrix into memory to calculate the statistics used for filtering.
2. The matrix is not continually held in memory in R workspace. So even if it is called into memory for a specific calculation, it will free up that memory once the calculation is done (assuming it can hold the full data in memory).



