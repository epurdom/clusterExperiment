#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// other_ids
// for a sample index i, determine what other indices in the same sample with it
// input:
// idx (integer) an index of a sample
// cl_vec (vector of integer) the clusterIds as returned by above in DList
// cl_len (vector of integer) the length of each cluster as returned by above in DList
// output:
//
// [[Rcpp::export]]
arma::uvec other_ids(arma::uword idx, arma::uvec cl_vec, arma::uvec cl_len) {



  arma::uvec ends = arma::cumsum(cl_len);
  // change below to something like: arma::uvec begins = arma::cumsum(c(1,head(clustLeng,-1)));
  arma::uvec begins = arma::cumsum(cl_len);

  arma::uword m = find(cl_vec == idx); //this may need to be a for loop?
  // TODO: check that sum(cl_vec == idx) == 1 -- if 0 return NA, if >1 error

  arma::uvec wh_cl = find(m <= ends && m >= begins);

  // check!
  // if(wh_cl.n_elem >1 || wh_cl.n_elem ==0) {
  //   stop("error in coding: finding range of clusterids")
  // }

  // retval = cl_vec[begins[wh_cl] : ends[wh_cl]];

  return(wh_cl);
}
