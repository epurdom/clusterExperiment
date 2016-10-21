# R package: clusterExperiment

Functions for running and comparing many different clusterings of single-cell sequencing data.

## Installation from bioconductor

We recommend installation of the package via bioconductor.

```r
source("https://bioconductor.org/biocLite.R")
biocLite("clusterExperiment")
```

To install the most recent version on the development branch of bioconductor, follow the above instructions, with the development version of bioconductor (see https://www.bioconductor.org/developers/how-to/useDevel/ for instructions [here](downloading development version of bioconductor) ).



## Installation of Github version:

While we generally try to keep the bioconductor devel version up-to-date with the master branch of the git repository, there is at times a lag between the two. You can install the github version via

```r
library(devtools)
install_github("epurdom/clusterExperiment")
```
## Status

| Resource:     |  Travis CI   |
| ------------- | ------------ |
| R CMD check master   | [![Build Status](https://travis-ci.org/epurdom/clusterExperiment.svg?branch=master)](https://travis-ci.org/epurdom/clusterExperiment)|
| ------------- | ------------ |
| R CMD check develop   | [![Build Status](https://travis-ci.org/epurdom/clusterExperiment.svg?branch=develop)](https://travis-ci.org/epurdom/clusterExperiment) |
| Test coverage |  [![Coverage Status](https://coveralls.io/repos/github/epurdom/clusterExperiment/badge.svg?branch=develop)](https://coveralls.io/github/epurdom/clusterExperiment?branch=develop) |

## Issues and bug reports

Please use https://github.com/epurdom/clusterExperiment/issues to submit issues, bug reports, and comments.
