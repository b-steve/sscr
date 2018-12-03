// Template to calculate the negative log-likelihood for a model with
// no random effects.
#include <TMB.hpp>
#include <fenv.h>
#include "detfns.h"
#include "utilities.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Reading in data.
  DATA_MATRIX(capt);
  DATA_IVECTOR(n_dets);
  DATA_MATRIX(mask_dists);
  DATA_INTEGER(n);
  DATA_INTEGER(n_traps);
  DATA_INTEGER(n_mask);
  DATA_SCALAR(mask_area);
  // Indicator for response type.
  DATA_INTEGER(resp_id);
  // Additional response parameters.
  DATA_VECTOR(resp_pars);
  // Indicator for detection function.
  DATA_INTEGER(detfn_id);
  // Indicators for parameter link functions.
  DATA_IVECTOR(link_ids);
  // Indicator for time-of-arrival data.
  DATA_INTEGER(toa_id);
  // Time-of-arrival data.
  DATA_MATRIX(toa_ssq);
  // Detection function parameter indices.
  DATA_IVECTOR(det_indices);
  // TOA parameter indices.
  DATA_INTEGER(toa_indices);
  // Declaring parameters.
  PARAMETER_VECTOR(link_pars);
  int n_pars = link_pars.size();
  // Back-transforming parameters.
  vector<Type> pars(n_pars);
  for (int i = 0; i < n_pars; i++){
    if (link_ids(i) == 0){
      pars(i) = exp(link_pars(i));
    } else if (link_ids(i) == 1){
      pars(i) = 1/(1 + exp(-link_pars(i)));
    }
  }
  ADREPORT(pars);
  // Extracting detection function parameters.
  int n_det_pars = det_indices.size();
  vector<Type> det_pars(n_det_pars);
  for (int i = 0; i < n_det_pars; i++){
    det_pars(i) = pars(det_indices(i));
  }
  // Extracting sigma_toa parameter.
  Type sigma_toa;
  if (toa_id == 1){
    sigma_toa = pars(toa_indices)/1000;
  }
  // Hazard rates for mask/trap combinations.
  matrix<Type> haz_mat(n_mask, n_traps);
  // Detection probabilities for mask/trap combinations.
  matrix<Type> prob_mat(n_mask, n_traps);
  // Detection probabilities for each mask point.
  vector<Type> prob_det(n_mask);
  // Generating hazard and probability matrices.
  for (int i = 0; i < n_mask; i++){
    for (int j = 0; j < n_traps; j++){
      prob_mat(i, j) = detfn(mask_dists(i, j), det_pars, detfn_id);
    }
  }
  haz_mat = prob_to_haz(prob_mat);
  // The sum of mask probabilities.
  Type sum_prob_det = 0;
  for (int i = 0; i < n_mask; i++){
    Type p_undet = Type(1);
    for (int j = 0; j < n_traps; j++){
      p_undet *= 1 - prob_mat(i, j);
    }
    // For binomial models, detection function is per *session*. Need
    // to account for this in overall detection probability.
    if (resp_id == 0){
      prob_det(i) = 1 - pow(p_undet, resp_pars(0));
    } else {
      prob_det(i) = 1 - p_undet;
    }
    sum_prob_det += prob_det(i);
  }
  // PMF for activity centres across the mask.
  vector<Type> f_loc(n_mask);
  f_loc = prob_det/sum_prob_det;
  // Likelihood contributions from capture histories.
  Type log_sum_integrands = 0;
  for (int i = 0; i < n; i++){
    Type integrand = 0;
    for (int j = 0; j < n_mask; j++){
      Type integrand_mask = 0;
      for (int k = 0; k < n_traps; k++){
	if (resp_id == 0){
	  integrand_mask += dbinom_sscr(capt(i, k), resp_pars(0), prob_mat(j, k), true);
	} else if (resp_id == 1){
	  integrand_mask += dpois_sscr(capt(i, k), haz_mat(j, k), true);
	}
      }
      // Time-of-arrival component.
      if (toa_id == 1){
	integrand_mask += (1 - n_dets(i))*log(sigma_toa) - (toa_ssq(i, j)/(2*pow(sigma_toa, 2)));
      }
      integrand += exp(integrand_mask);
    }
    log_sum_integrands += log(integrand + DBL_MIN);
  }
  Type f = -log_sum_integrands;
  Type esa = mask_area*sum_prob_det;
  ADREPORT(esa);
  Type D = n/esa;
  ADREPORT(D);
  // Extra bit that falls out of log-likelihood.
  f -= -n*log(sum_prob_det);
  // Contribution from number of animals detected.
  f -= dpois_sscr(n, D*esa, true);
  return f;
}
