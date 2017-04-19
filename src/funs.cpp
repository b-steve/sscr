#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// Stolen from the ascr package.
// [[Rcpp::export]]
NumericMatrix make_toa_ssq(const NumericMatrix& capt, const NumericMatrix& dists, const double& sound_speed){
  int n = capt.nrow();
  int n_traps = capt.ncol();
  int n_mask = dists.ncol();
  NumericMatrix out(n, n_mask);
  int n_dets;
  int index;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n_mask; j++){
      n_dets = 0;
      for (int k = 0; k < n_traps; k++){
	if (capt(i, k) > 0) n_dets++;
      }
      NumericVector delts(n_dets);
      index = 0;
      for (int k = 0; k < n_traps; k++){
	if (capt(i, k) > 0){
	  delts(index) = capt(i, k) - dists(k, j)/sound_speed;
	  index++;
	}
	out(i, j) = sum(pow(delts - mean(delts), 2));
      }
    }
  }  
  return out;
}
