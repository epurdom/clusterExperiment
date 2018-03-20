// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <Rcpp.h>
using namespace Rcpp;
// Search pairs of samples that co-cluster across subsamples
//
// Assume that our input is a matrix, with N columns and B rows
// (the number of subsamples), storing integers -- the cluster labels.
//
// The output is a matrix with the co-clusters, but only the lower triangle
// is populated.
// [[Rcpp::export]]
NumericMatrix search_pairs(IntegerMatrix clusterings) {

  // N columns; B rows
  int N = clusterings.ncol();
  NumericMatrix retval(N, N);

  for(int i = 1; i < clusterings.ncol(); i++) {

    for(int j = 0; j < i; j++) {

      double s = sum(na_omit(clusterings.column(i) == clusterings.column(j)));

      double tot_na = sum(!(is_na(clusterings.column(i)) | is_na(clusterings.column(j))));

      retval(i, j) = s / tot_na;

    }

  }

  return retval;
}

