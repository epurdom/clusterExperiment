# R package: clusterExperiment

Functions for running and comparing many different clusterings of single-cell sequencing data.

## News and Updates

* Note that some important bugs were fixed since the Bioconductor release. They have been updated in the release version of Bioconductor (version 2.0.2) as well as the development version. 

* Version 2.0.0 is now available on Bioconductor with many new changes. Checkout out a [brief description](https://github.com/epurdom/clusterExperiment/blob/master/update2.0.md) of the major changes. For complete details, see the [NEWS](https://github.com/epurdom/clusterExperiment/blob/master/NEWS) file.

  The new 2.0.0 version does make use of packages only available in the new release of bioconductor (and thus the new release of R).

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
