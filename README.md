# R package: clusterExperiment

Functions for running and comparing many different clusterings of single-cell sequencing data.

## Installation From Bioconductor

We recommend installation of the package via bioconductor.

```r
source("https://bioconductor.org/biocLite.R")
biocLite("clusterExperiment")
```

To install the most recent version on the development branch of bioconductor, follow the above instructions, with the development version of bioconductor (see  [here](https://www.bioconductor.org/developers/how-to/useDevel/) for instructions).

## Installation of Github Version:

We generally try to keep the bioconductor *devel* version up-to-date with the *master* branch of this git repository, but there can be at times a lag between the two. You can install the github version via

```r
library(devtools)
install_github("epurdom/clusterExperiment")
```

## Development branch:

The `develop` branch is our development branch where we are actively updating features, and may contain bugs or be in the process of being updated. You should not use the `develop` branch unless it passes TravisCI checks (see below) and you want to be using a *very* beta version.

The development branch can be installed via the `install_github` command above, but indicating the `develop` branch:

```r
library(devtools)
install_github("epurdom/clusterExperiment", ref="develop")
```

## Status

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
