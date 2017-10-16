# R package: clusterExperiment

Functions for running and comparing many different clusterings of single-cell sequencing data.

## Installation from bioconductor

We recommend installation of the package via bioconductor.

```r
source("https://bioconductor.org/biocLite.R")
biocLite("clusterExperiment")
```

To install the most recent version on the development branch of bioconductor, follow the above instructions, with the development version of bioconductor (see  [here](https://www.bioconductor.org/developers/how-to/useDevel/) for instructions).

## Installation of Github version:

We generally try to keep the bioconductor devel version up-to-date with the master branch of the git repository, but there can be at times a lag between the two. You can install the github version via

```r
library(devtools)
install_github("epurdom/clusterExperiment")
```

## Development branch:

The `develop` branch is our development branch where we are actively updating features, and may contain bugs. You should not use the `develop` branch unless it passes TravisCI checks and you want to be using a *very* beta version.

## Status

Below are status checks for the package. Note that occassionally errors do not appear here immediately. Clicking on the link will give you the most up-to-date status.

| Resource:     |  Travis CI   |
| ------------- | ------------ |
| R CMD check master   | [![Build Status](https://travis-ci.org/epurdom/clusterExperiment.svg?branch=master)](https://travis-ci.org/epurdom/clusterExperiment)|
| ------------- | ------------ |
| R CMD check develop   | [![Build Status](https://travis-ci.org/epurdom/clusterExperiment.svg?branch=develop)](https://travis-ci.org/epurdom/clusterExperiment) |
| Test coverage |  [![Coverage Status](https://coveralls.io/repos/github/epurdom/clusterExperiment/badge.svg?branch=develop)](https://coveralls.io/github/epurdom/clusterExperiment?branch=develop) |
| ------------- | ------------ |

| Resource:     |  Bioconductor  |
| ------------- | ------------ |
| Bioc Release  | [![BiocDevel Status](http://bioconductor.org/shields/build/release/bioc/clusterExperiment.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/clusterExperiment/)|
| ------------- | ------------ |
| Bioc Development  | [![BiocDevel Status](http://bioconductor.org/shields/build/devel/bioc/clusterExperiment.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/clusterExperiment/)|
| ------------- | ------------ |

## Issues and bug reports

Please use https://github.com/epurdom/clusterExperiment/issues to submit issues, bug reports, and comments.
