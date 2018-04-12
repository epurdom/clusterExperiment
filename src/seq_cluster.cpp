// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <algorithm>
#include <iostream>
using namespace Rcpp;

// Function to calculate the stability across the given combination
//	y is a combination (row of index_m) giving clusters to compare stability from k, k+1
float calc_beta(IntegerVector y, List candidates, int seq_num, String beta_num) {
  // written generally enough to deal with seq_num>2; could be a lot simpler with seq.num=2.

  List temp(seq_num);

  for(int z=0; z<seq_num; z++) {
    List each_can = candidates[z];
    temp[z] = each_can[y[z]];
  }

  IntegerVector itemp = temp[0];
  IntegerVector utemp;

  if(beta_num == "all" || beta_num == "first") {
    utemp = temp[0];
  }

  if(beta_num == "last") {
    utemp = temp[seq_num - 1];
  }

  for(int j=1; j<seq_num; j++) {
    IntegerVector jtemp = temp[j];
    itemp = intersect(itemp, jtemp);

    if(beta_num == "all") {
      utemp = union_(utemp, jtemp);
    }
  }

  float retval =  (float)itemp.size() / (float)utemp.size();

  return(retval);

}



//' Sequential clustering
//'
//' Given a data matrix, this function will call clustering
//' routines, and sequentially remove best clusters, and iterate to find
//' clusters.
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List do_seq_cluster(arma::mat& x, arma::mat& diss, int N, int k0,
                    int remain_n, int seq_num, int k_min, int k_max,
                    int top_can, float beta,
                    IntegerMatrix index_m, String beta_num,
                    int wh_return, String input,
                    Function updateClustering) {

  // setup
  int remain = N; //keep track of how many samples not yet clustered (stop when less than remain_n)
  int nfound = 0; //keep track of how many clusters found/removed so far
  bool found = true; //has a cluster been found/removed in last iteration
  int k_start = k0; //the starting k for the next cluster
  int k = k0;
  int why_stop=0;

  List candidates(seq_num); //list of length seq_num of possible clusters found for each k to be compared
  List tclust; //list of final cluster identifications (indices of rows of x)
  IntegerVector kstart; //the starting k for the cluster
  IntegerVector kend; //the ending k for the cluster
  int current_start;
  int iter=1;
  while(remain >= remain_n && (found || k <= k_max)) {
    std::cout << "itearation: " << iter << std::endl;
    if(found) {
      k = k_start;
      current_start = k_start;
      for(int i=0; i<seq_num; i++) {
        int newk = k + i;
        List res = updateClustering(newk, x, diss);
        int s = res.size();
        if(s > 0) {
          IntegerVector idx = seq(0, std::min(top_can, s) - 1);
          res = res[idx];
        }
        candidates[i] = res;
      }
    } else {
      IntegerVector idx = seq_len(candidates.size());
      candidates = candidates[idx != 1];

      int newk = k + seq_num - 1;
      List res = updateClustering(newk, x, diss);
      int s = res.size();
      if(s > 0) {
        IntegerVector idx = seq(0, std::min(top_can, s) - 1);
        res = res[idx];
      }

      candidates.push_back(res);

    }

    // check whether all got top_can values for each -- could be less.
    int nClusterPerK;
    LogicalMatrix is_invalid(index_m.nrow(), seq_num);

    for(int i=0; i<seq_num; i++) {
      List temp = candidates[i];
      nClusterPerK = temp.size();
      is_invalid.column(i) = index_m.column(i) > (nClusterPerK - 1);
    }

    int ninv=0;
    for(int i=0; i<index_m.nrow(); i++) {
      ninv = ninv + is_true(any(is_invalid.row(i)));
    }

    IntegerMatrix index_red(index_m.nrow() - ninv, seq_num);

    // all invalid -- probably means that for some k there were no candidates found.
    // So should stop.
    if(ninv == index_m.nrow()) {
      why_stop = 1;
      break;
    }

    if(ninv > 0) {
      int j=0;
      for(int i=0; i<index_m.nrow(); i++) {
        if(is_true(any(is_invalid.row(i)))) {
          continue;
        }
        index_red.row(j) = index_m.row(i);
        j++;
      }
    } else {
      index_red = index_m;
    }

    // Calculate the stability pairwise between all of cluster combinations
    NumericVector beta_temp(index_red.nrow());
    for(int i=0; i<index_red.nrow(); i++) {
      beta_temp[i] = calc_beta(index_red.row(i), candidates, seq_num, beta_num);
    }

    std::cout << "beta_temp:" << beta_temp << std::endl;

    if(is_true(any(beta_temp >= beta))) {
      found = true;

      if(k_start > k_min) {
        k_start = k_start - 1;
      }

      List list_temp = candidates[wh_return];
      IntegerVector found_temp = list_temp[index_red(which_max(beta_temp), wh_return)];

      std::cout << "Found cluster number " << nfound + 1 << std::endl;
      std::cout << "Cluster size: " << found_temp.size()  << std::endl;

      kend.push_back(k + seq_num - 1);
      kstart.push_back(current_start);
      tclust.push_back(found_temp); // the indexes of found_temp are wrong; need sample names!
      nfound++;

      IntegerVector all_idx = seq(1, x.n_cols);
      IntegerVector tofind = setdiff(all_idx, found_temp) - 1;
      arma::uvec to_find = as<arma::uvec>(tofind);

      if(input == "X") {
        x = x.cols(to_find);
        diss = diss.cols(to_find);
        diss = diss.rows(to_find);
      }

      if(input == "diss") {
        x = x.cols(to_find);
        diss = diss.cols(to_find);
        diss = diss.rows(to_find);
      }

      remain = remain - found_temp.size();
      std::cout << "Samples remaining: " << remain << std::endl;

    } else {
      found = false;
      k++;
    }

    iter++;
  }

  return List::create(Named("tclust") = tclust,
                      Named("nfound") = nfound,
                      Named("kstart") = kstart,
                      Named("kend") = kend,
                      Named("remain") = remain,
                      Named("remain.n") = remain_n,
                      Named("found") = found,
                      Named("k") = k,
                      Named("k.max") = k_max,
                      Named("whyStop") = why_stop);

}
