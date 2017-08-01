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
  // Capture histories.
  DATA_MATRIX(capt);
  // Number of detections.
  DATA_IVECTOR(n_dets);
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
  // Indicator for detection function.
  DATA_INTEGER(detfn_id);
  // Indicator for detection function scale.
  DATA_INTEGER(detfn_scale_id);
  // Indicator for dependence structure.
  DATA_INTEGER(cov_id);
  // Indicator for random effects scale.
  DATA_INTEGER(re_scale_id);
  // Detection probabilities (from cov_detprob.cpp).
  DATA_VECTOR(det_probs);
  // Indicator for time-of-arrival data.
  DATA_INTEGER(toa_id);
  // Time-of-arrival data.
  DATA_MATRIX(toa_ssq);
  // Detection function parameters.
  DATA_VECTOR(det_pars);
  // Covariance parameters.
  DATA_VECTOR(cov_pars);
  // Time-of-arrival paramter.
  DATA_SCALAR(sigma_toa);
  // Indicator for conditional likelihood.
  DATA_INTEGER(conditional_n);
  // Density parameter. Required if likelihood is not conditional.
  DATA_SCALAR(D);
  // Latent variables.
  PARAMETER_MATRIX(u);
  // Hazard rates for mask/trap combinations.
  matrix<Type> haz_mat(n_mask, n_traps);
  // Detection probabilities for mask/trap combinations.
  matrix<Type> prob_mat(n_mask, n_traps);
  if (detfn_scale_id == 0){
    for (int i = 0; i < n_mask; i++){
      for (int j = 0; j < n_traps; j++){
	haz_mat(i, j) = detfn(mask_dists(i, j), det_pars, detfn_id);
      }
    }
    prob_mat = haz_to_prob(haz_mat);
  } else if (detfn_scale_id == 1){
    for (int i = 0; i < n_mask; i++){
      for (int j = 0; j < n_traps; j++){
	prob_mat(i, j) = detfn(mask_dists(i, j), det_pars, detfn_id);
      }
    }
    haz_mat = prob_to_haz(prob_mat);
  }
  // The sum of mask probabilities.
  Type sum_det_probs = 0;
  for (int i = 0; i < n_mask; i++){
    sum_det_probs += det_probs(i);
  }
  Type u_use;
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
      Type integrand_mask = 0;
      for (int k = 0; k < n_traps; k++){
	if (cov_id == 3){
	  u_use = u(i, 0);
	} else {
	  u_use = u(i, k);
	}
	// Expected counts and probabilities.
	Type e_count = exp(log(haz_mat(j, k) + 1e-12) + u_use) + DBL_MIN;
	Type e_prob = 1 - exp(-e_count);
	if (re_scale_id == 0){
	  e_count = exp(log(haz_mat(j, k) + 1e-12) + u_use) + DBL_MIN;
	  e_prob = haz_to_prob(e_count);
	} else if (re_scale_id == 1){
	  e_prob = invlogit(logit(prob_mat(j, k) + 1e-12) + u_use) + DBL_MIN;
	  e_count = prob_to_haz(e_prob);
	}
	if (resp_id == 0){
	  integrand_mask += dbinom_sscr(capt(i, k), resp_pars(0), e_prob, true);
	} else if (resp_id == 1){
	  integrand_mask += dpois_sscr(capt(i, k), e_count, true);
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
  f -= log_sum_integrands;
  // Extra bit that falls out of log-likelihood.
  f -= -n*log(sum_det_probs);
  // Likelihood component due to n.
  if (conditional_n == 0){
    std::cout << "f before: " << f << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "pmf(n): " << dpois_sscr(Type(n), Type(D*mask_area*sum_det_probs), true) << std::endl;
    f -= dpois_sscr(Type(n), Type(D*mask_area*sum_det_probs), true);
    std::cout << "f after: " << f << std::endl;
  }
  if (cov_id == 3){
    for (int i = 0; i < n; i++){
      f -= dnorm(u(i, 0), Type(0), cov_pars(0), true);
    }
  } else {
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
	      // Total dependence (individual-level random effects).
	      sigma_u_mat(j, k) = pow(cov_pars(0), 2);
	      sigma_u_mat(k, j) = pow(cov_pars(0), 2);
	    } else if (cov_id == 4){
	      // Linear combination of exponential covariance functions.
	      sigma_u_mat(j, k) = pow(cov_pars(0), 2)*(cov_pars(1)*exp(-cov_pars(2)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(2)*trap_dists(j, k)), 2) + (1 - cov_pars(1))*exp(-cov_pars(3)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(3)*trap_dists(j, k)), 2));
	      sigma_u_mat(k, j) = pow(cov_pars(0), 2)*(cov_pars(1)*exp(-cov_pars(2)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(2)*trap_dists(j, k)), 2) + (1 - cov_pars(1))*exp(-cov_pars(3)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(3)*trap_dists(j, k)), 2));
	    } else if (cov_id == 5){
	      // Squared exponential covariance function.
	      sigma_u_mat(j, k) = pow(cov_pars(0), 2)*exp(-pow(trap_dists(j, k), 2)/pow(cov_pars(1), 2));
	      sigma_u_mat(k, j) = pow(cov_pars(0), 2)*exp(-pow(trap_dists(j, k), 2)/pow(cov_pars(1), 2));
	    } else {
	      exit(1111);
	    }
	  }
	}
      }
      // Contribution from latent variables (note MVNORM returns the
      // negative-log of the density).
      f += MVNORM(sigma_u_mat)(u.row(i));
    }
  }
  return f;
}
