# R package: clusterExperiment

Functions for running and comparing many different clusterings of single-cell sequencing data.

## News and Updates

Version 2.0.0 is now available on with many new changes. Checkout out a [brief description](https://github.com/epurdom/clusterExperiment/blob/master/update2.0.md) of the major changes. For complete details, see the [NEWS](https://github.com/epurdom/clusterExperiment/blob/master/NEWS) file.

The new 2.0.0 version does make use of packages only available in the new release of bioconductor (and thus the new release of R).

## Publications

The paper acompanying this package can be found at:

Davide Risso, Liam Purvis, Russell Fletcher, Diya Das, John Ngai, Sandrine Dudoit, and Elizabeth Purdom (March 2018) "clusterExperiment and RSEC: A Bioconductor package and framework for clustering of single-cell and other large gene expression datasets" bioRxiv, (https://www.biorxiv.org/content/early/2018/03/12/280545)[https://www.biorxiv.org/content/early/2018/03/12/280545].

There is a github repository ((epurdom/RSECPaper)[https://github.com/epurdom/RSECPaper] for this paper that gives the code for reproducing the analysis in that manuscript.

A compiled version of the vignette (i.e. tutorial) of the github version of `clusterExperiment` (corresponding to the `master` branch) can be found (here)[http://bioconductor.org/packages/devel/bioc/vignettes/clusterExperiment/inst/doc/clusterExperimentTutorial.html] 

The compiled version of the version of `clusterExperiment` available with the latest release of Bioconductor can be found (here)[http://bioconductor.org/packages/release/bioc/vignettes/clusterExperiment/inst/doc/clusterExperimentTutorial.html] 

## Installation instructions

### Installation From Bioconductor

We recommend installation of the package via bioconductor.

```r
source("https://bioconductor.org/biocLite.R")
biocLite("clusterExperiment")
```

To install the most recent version on the development branch of bioconductor, follow the above instructions, with the development version of bioconductor (see  [here](https://www.bioconductor.org/developers/how-to/useDevel/) for instructions).

### Installation of Github Version:

We generally try to keep the bioconductor *devel* version up-to-date with the *master* branch of this git repository. But since this can require installing R development version, it can be convenient to be able to get just the most recent version from github. 

You can install the github version via

```r
library(devtools)
install_github("epurdom/clusterExperiment")
```

### Development branch:

The `develop` branch is our development branch where we are actively updating features, and may contain bugs or be in the process of being updated. You should not use the `develop` branch unless it passes TravisCI checks (see below) and you want to be using a *very* beta version.

The development branch can be installed via the `install_github` command above, but indicating the `develop` branch:

```r
library(devtools)
install_github("epurdom/clusterExperiment", ref="develop")
```

## Status Checks

Below are status checks for the package. Note that occassionally errors do not appear here immediately. Clicking on the link will give you the most up-to-date status.

| Resource:     |  Status   |
| ------------- | ------------ |
| Bioc Release  | [![BiocDevel Status](http://bioconductor.org/shields/build/release/bioc/clusterExperiment.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/clusterExperiment/)|
| Bioc Development  | [![BiocDevel Status](http://bioconductor.org/shields/build/devel/bioc/clusterExperiment.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/clusterExperiment/)|
| Travis CI master   | [![Build Status](https://travis-ci.org/epurdom/clusterExperiment.svg?branch=master)](https://travis-ci.org/epurdom/clusterExperiment) |
| Travis CI develop   | [![Build Status](https://travis-ci.org/epurdom/clusterExperiment.svg?branch=develop)](https://travis-ci.org/epurdom/clusterExperiment) |
| Appveyor master | [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/epurdom/clusterExperiment?branch=master&svg=true)](https://ci.appveyor.com/project/epurdom/clusterExperiment) |
| Appveyor develop | [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/epurdom/clusterExperiment?branch=develop&svg=true)](https://ci.appveyor.com/project/epurdom/clusterExperiment) |
| Test coverage |  [![Coverage Status](https://coveralls.io/repos/github/epurdom/clusterExperiment/badge.svg?branch=develop)](https://coveralls.io/github/epurdom/clusterExperiment?branch=develop) |

## Issues and bug reports

Please use https://github.com/epurdom/clusterExperiment/issues to submit issues, bug reports, and comments.
