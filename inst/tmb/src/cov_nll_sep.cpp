// Template to calculate the negative log-likelihood for a model with
// independent random effects.
#include <TMB.hpp>
#include <fenv.h>
#include "utilities.h"
#include "detfns.h"
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Capture histories.
  DATA_VECTOR(capt);
  // Number of detections.
  DATA_INTEGER(n_dets);
  // Distance from mask point to all traps.
  DATA_MATRIX(mask_dists);
  // Distances between traps.
  DATA_MATRIX(trap_dists);
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
  // Indicator for dependence structure.
  DATA_INTEGER(cov_id);
  // Indicator for random effects scale.
  DATA_INTEGER(re_scale_id);
  // Detection probabilities (from cov_detprob.cpp).
  DATA_VECTOR(det_probs);
  // Indicator for time-of-arrival data.
  DATA_INTEGER(toa_id);
  // Time-of-arrival data.
  DATA_VECTOR(toa_ssq);
  // Indicators for detection parameter link functions.
  DATA_IVECTOR(link_det_ids);
  // Number of detection function parameters.
  int n_det_pars = link_det_ids.size();
  // Indicators for covariance parameter link functions.
  DATA_IVECTOR(link_cov_ids);
  // Number of covariance function parameters.
  int n_cov_pars = link_cov_ids.size();
  // Detection function parmaters.
  PARAMETER_VECTOR(link_det_pars);
  // Covariance parameters.
  PARAMETER_VECTOR(link_cov_pars);
  // Time-of-arrival parameter.
  PARAMETER(link_sigma_toa);
  // Density parameter.
  PARAMETER(link_D);
  // Latent variables.
  PARAMETER_VECTOR(u);
  // Back-transforming detection function parameters.
  vector<Type> det_pars(n_det_pars);
  for (int i = 0; i < n_det_pars; i++){
    if (link_det_ids(i) == 0){
      det_pars(i) = exp(link_det_pars(i));
    } else if (link_det_ids(i) == 1){
      det_pars(i) = 1/(1 + exp(-link_det_pars(i)));
    }
  }
  // Back-transforming covariance function parameters.
  vector<Type> cov_pars(n_cov_pars);
  for (int i = 0; i < n_cov_pars; i++){
    if (link_cov_ids(i) == 0){
      cov_pars(i) = exp(link_cov_pars(i));
    } else if (link_cov_ids(i) == 1){
      cov_pars(i) = 1/(1 + exp(-link_cov_pars(i)));
    }
  }
  // Back-transforming sigma_toa, including readjustment to ms.
  Type sigma_toa = exp(link_sigma_toa)/1000;
  // Hazard rates for mask/trap combinations.
  matrix<Type> haz_mat(n_mask, n_traps);
  // Detection probabilities for mask/trap combinations.
  matrix<Type> prob_mat(n_mask, n_traps);
  for (int i = 0; i < n_mask; i++){
    for (int j = 0; j < n_traps; j++){
      prob_mat(i, j) = detfn(mask_dists(i, j), det_pars, detfn_id);
    }
  }
  haz_mat = prob_to_haz(prob_mat);
  // The sum of mask probabilities.
  Type sum_det_probs = 0;
  for (int i = 0; i < n_mask; i++){
    sum_det_probs += det_probs(i);
  }
  // Declaring the vector of random effects to use.
  Type u_use;
  // Likelihood contributions from capture histories.
  Type integrand = 0;
  for (int j = 0; j < n_mask; j++){
    Type integrand_mask = 0;
    for (int k = 0; k < n_traps; k++){
      if (cov_id == 6){
	u_use = 0;
      } else if (cov_id == 3){
	u_use = u(0);
      } else {
	u_use = u(k);
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
	integrand_mask += dbinom_sscr(capt(k), resp_pars(0), e_prob, true);
      } else if (resp_id == 1){
	integrand_mask += dpois_sscr(capt(k), e_count, true);
      }
    }
    // Time-of-arrival component.
    if (toa_id == 1){
      integrand_mask += (1 - n_dets)*log(sigma_toa) - (toa_ssq(j)/(2*pow(sigma_toa, 2)));
    }
    integrand += exp(integrand_mask);
  }
  Type f = -log(integrand + DBL_MIN);
  // Extra bit that falls out of log-likelihood.
  f -= -log(sum_det_probs);
  // Sorting out spatial random effects.
  if (cov_id == 3){
      f -= dnorm(u(0), Type(0), cov_pars(0), true);
  } else if (cov_id != 6){
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
	    // NOTE: This is actually handled above.
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
    f += MVNORM(sigma_u_mat)(u);
  }
  return f;
}