// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <Rcpp.h>
using namespace Rcpp;
//' Search pairs of samples that co-cluster across subsamples
//'
// Assume that our input is a matrix, with N rows and B columns
// (the number of subsamples), storing integers -- the cluster labels.
//
// The output should be a list of length N - 1 with the proportion of times that
// an observation is co-clustered with the ones before it.
// Note that NA are translated to 0 if umat -- need to work with mat.
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





// arma::umat search_pairs(arma::mat clusterings) {
//
//   // For each column, store the number of 0's (NA's)
//   int N = clusterings.n_cols;
//   // arma::uvec tot_obs(N);
//   // for(int i = 1; i < clusterings.n_cols; i++) {
//   //
//   //   tot_obs[i] = sum(clusterings.col(i).is_finite);
//   //
//   // }
//   //
//   // std::cout << (tot_obs) << std::endl;
//
//   // N columns; B rows
//   arma::umat retval(N, N, arma::fill::zeros);
//
//   for(int i = 1; i < clusterings.n_cols; i++) {
//
//     for(int j = 0; j < i; j++) {
//
//       arma::uword s = sum(clusterings.col(i) == clusterings.col(j));
//
//       // arma::uword tot = arna::min(tot_i, tot_j)
//       // int num_obs = clusterings.n_rows - arma::max(arma::sum(arma::uvec::is_na(clusterings.col(i))),
//       //                                             arma::sum(arma::uvec::is_na(clusterings.col(j))))
//       // double prop = s / num_obs;
//
//       LogicalVector = IntegerVector.is_na(clusterings.col(i))
//
//         std::cout << (clusterings.col(i) == clusterings.col(j)) << std::endl;
//       // std::cout << (finite_i) << std::endl;
//       // std::cout << (finite_j) << std::endl;
//
//       //std::cout << (s) << std::endl;
//
//       retval(i, j) = s;
//
//     }
//
//   }
//
//   return retval;
// }

// IntegerMatrix search_pairs(IntegerMatrix clusterings) {
//
//   // N columns; B rows
//   int N = clusterings.ncol();
//   IntegerMatrix retval(N, N);
//
//   for(int i = 1; i < clusterings.ncol(); i++) {
//
//     for(int j = 0; j < i; j++) {
//
//       LogicalVector idx = clusterings.column(i) == clusterings.column(j);
//       LogicalVector na_idx = IntegerV ector::is_na(clusterings.column(i));
//
//       int s = sum(idx);
//       // int num_obs = clusterings.n_rows - arma::max(arma::sum(arma::uvec::is_na(clusterings.col(i))),
//       //                                             arma::sum(arma::uvec::is_na(clusterings.col(j))))
//       // double prop = s / num_obs;
//
//       std::cout << (na_idx) << std::endl;
//
//       retval(i, j) = s;
//
//     }
//
//   }
//
//   return retval;
// }


// int search_pairs(arma::umat clusterings) {
//
//   // loop for each sample and observation
//   for(arma::uword cl=0; cl < clusterings.n_cols; ++cl) {
//
//
//     clusterIds <- unname(unlist(tapply(1:N, classX, identity, simplify=FALSE)))
//
//     whHave<-which(sapply(clusterList,function(ll){ii%in%ll$clusterIds}))
//
//     auto sample_end = clusterings.end_col(cl);
//
//     for(auto sample = clusterings.begin_col(cl); sample != sample_end; ++sample) {
//
//       std::cout << (*sample) << std::endl;
//
//       // get list of those indices sample it was sampled
//       clusterIds =
//       // calculate number of times sampled with (denominator)
//
//       // get those indices clustered with and tabulate
//
//       // return proportions of times it is together with the others
//
//
//     }
//
//   }
//
//
//   for (int i = 0; i < clusterings.n_rows; i++){
//     for (int j = 0; j < clusterings.n_cols; j++){
//       int one = clusterings(i,j);
//
//       std::cout << (one) << std::endl;
//
//     }
//   }
//
//   return 0;
// }















// other_ids
// for a sample index i, determine what other indices in the same sample with it
// input:
// idx (integer) an index of a sample
// cl_vec (vector of integer) the clusterIds as returned by above in DList
// cl_len (vector of integer) the length of each cluster as returned by above in DList
// output:
//
// arma::uvec other_ids(arma::uword idx, arma::uvec cl_vec, arma::uvec cl_len) {
//
//   arma::uvec m = find(cl_vec == idx);
//   // TODO: check that sum(cl_vec == idx) == 1 -- if 0 return NA, if >1 error
//
//   arma::uvec ends = arma::cumsum(cl_len);
//   // change below to something like: arma::uvec begins = arma::cumsum(c(1,head(clustLeng,-1)));
//   arma::uvec begins = arma::cumsum(cl_len);
//
//
//   arma::uvec wh_cl = find(m[0] <= ends && m[0] >= begins);
//
//   // check!
//   // if(wh_cl.n_elem >1 || wh_cl.n_elem ==0) {
//   //   stop("error in coding: finding range of clusterids")
//   // }
//
//   // retval = cl_vec[begins[wh_cl] : ends[wh_cl]];
//
//   return(wh_cl);
// }
