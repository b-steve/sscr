// Template to calculate the negative log-likelihood for a model with
// independent random effects.
#include <TMB.hpp>
#include <fenv.h>
#include "detfns.h"
#include "utilities.h"
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Capture histories
  DATA_MATRIX(capt);
  // Distance from mask point to all traps.
  DATA_MATRIX(mask_dists);
  // Distances between traps.
  DATA_MATRIX(trap_dists);
  // Number of detected indivduals.
  DATA_INTEGER(n);
  // Number of traps.
  DATA_INTEGER(n_traps);
  // Number of mask points.
  DATA_INTEGER(n_mask);
  // Area of each mask point.
  DATA_SCALAR(mask_area);
  // Indicator for response distribution.
  DATA_INTEGER(resp_id);
  // Additional response parameters.
  DATA_VECTOR(resp_pars);
  // Indicator for detection function ID.
  DATA_INTEGER(detfn_id);
  // Indicator for dependence structure.
  DATA_INTEGER(cov_id);
  // Detection probabilities (from cov_detprob.cpp).
  DATA_VECTOR(det_probs);
  // Detection function parameters.
  DATA_VECTOR(det_pars);
  // Covariance parameters.
  DATA_VECTOR(cov_pars);
  // Latent variables.
  PARAMETER_MATRIX(u);
  // Hazard rates for mask/trap combinations.
  matrix<Type> haz_mat(n_mask, n_traps);
  // The sum of mask probabilities.
  Type sum_det_probs = 0;
  for (int i = 0; i < n_mask; i++){
    for (int j = 0; j < n_traps; j++){
      // Calculating encounter rate.
      haz_mat(i, j) = detfn(mask_dists(i, j), det_pars, detfn_id);
    }
    sum_det_probs += det_probs(i);
  }
  // PMF for activity centres across the mask.
  vector<Type> f_loc(n_mask);
  f_loc = det_probs/sum_det_probs;
  // Joint density of data and latent variables.
  Type f = 0;
  // Likelihood contributions from capture histories.
  Type log_sum_integrands = 0;
  for (int i = 0; i < n; i++){
    Type integrand = 0;
    for (int j = 0; j < n_mask; j++){
      Type integrand_mask = 1;
      for (int k = 0; k < n_traps; k++){
	Type e_count = exp(log(haz_mat(j, k)) + u(i, k)) + DBL_MIN;
	if (resp_id == 0){
	  Type e_prob = 1 - exp(-e_count);
	  integrand_mask *= dbinom_sscr(capt(i, k), resp_pars(0), e_prob, false);
	} else if (resp_id == 1){
	  integrand_mask *= dpois(capt(i, k), e_count, false);
	}
      }
      integrand += integrand_mask;
    }
    log_sum_integrands += log(integrand + DBL_MIN);
  }
  f -= log_sum_integrands;
  // Extra bit that falls out of log-likelihood.
  f -= -n*log(sum_det_probs);
  for (int i = 0; i < n; i++){
    // Variance-covariance matrix for latent variables.
    matrix<Type> sigma_u_mat(n_traps, n_traps);
    for (int j = 0; j < n_traps; j++){
      for (int k = j; k < n_traps; k++){ 
	if (j == k){
	  sigma_u_mat(j, k) = pow(cov_pars(0), 2);
	} else {
	  if (cov_id == 0){
	    // Independent random effects
	    sigma_u_mat(j, k) = 0;
	    sigma_u_mat(k, j) = 0;
	  } else if (cov_id == 1){
	    // Exponential covariance function.
	    sigma_u_mat(j, k) = pow(cov_pars(0), 2)*exp(-trap_dists(j, k)/cov_pars(1));
	    sigma_u_mat(k, j) = pow(cov_pars(0), 2)*exp(-trap_dists(j, k)/cov_pars(1));
	  } else if (cov_id == 2){
	    // Matern covariance function.
	  } else if (cov_id == 3){
	    // Total dependence.
	    sigma_u_mat(j, k) = pow(cov_pars(0), 2);
	    sigma_u_mat(k, j) = pow(cov_pars(0), 2);
	  } else if (cov_id == 4){
	    sigma_u_mat(j, k) = pow(cov_pars(0), 2)*(cov_pars(1)*exp(-cov_pars(2)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(2)*trap_dists(j, k)), 2) + (1 - cov_pars(1))*exp(-cov_pars(3)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(3)*trap_dists(j, k)), 2));
	    sigma_u_mat(k, j) = pow(cov_pars(0), 2)*(cov_pars(1)*exp(-cov_pars(2)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(2)*trap_dists(j, k)), 2) + (1 - cov_pars(1))*exp(-cov_pars(3)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(3)*trap_dists(j, k)), 2));
	  }
	}
      }
    }
    // Contribution from latent variables (note MVNORM returns the
    // negative-log of the density).
    f += MVNORM(sigma_u_mat)(u.row(i));
  }
  return f;
}
