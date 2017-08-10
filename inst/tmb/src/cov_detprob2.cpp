// Template to calculate detection probability at a particular mask
// point for a model with some covariance function for the random
// effects.
#include <TMB.hpp>
#include <fenv.h>
#include "detfns.h"
#include "utilities.h"
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Distance from mask point to all traps.
  DATA_VECTOR(mask_dists);
  // Distances between traps.
  DATA_MATRIX(trap_dists);
  // Number of traps.
  DATA_INTEGER(n_traps);
  // Indicator for detection function.
  DATA_INTEGER(detfn_id);
  // Indicator for detection function scale.
  DATA_INTEGER(detfn_scale_id);
  // Indicator for response type.
  DATA_INTEGER(resp_id);
  // Additional response parameters.
  DATA_VECTOR(resp_pars)
  // Indicator for dependence structure.
  DATA_INTEGER(cov_id);
  // Indicator for random effects scale.
  DATA_INTEGER(re_scale_id);
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
  // Back-transforming sigma_toa.
  Type sigma_toa = exp(link_sigma_toa);
  // Back-transforming density parameter.
  Type D = exp(link_D);
  // Declaring latent variables.
  Type u_use;
  // Overall probability of nondetection.
  Type p_total_evade = 1;
  // Negative-log joint density of probability of detection and latent variables.
  Type f = 0;
  // Baseline hazards and probabilities.
  Type base_haz;
  Type base_prob;
  // Actual hazards and probabilities.
  Type haz;
  Type prob;
  // Probability of capture.
  for (int i = 0; i < n_traps; i++){
    // Calculating baseline hazard rate and probability.
    if (detfn_scale_id == 0){
      base_haz = detfn(mask_dists(i), det_pars, detfn_id);
      base_prob = haz_to_prob(base_haz);
    } else if (detfn_scale_id == 1){
      base_prob = detfn(mask_dists(i), det_pars, detfn_id);
      base_haz = prob_to_haz(base_prob);
    }
    if (cov_id == 3){
      u_use = u(0);
    } else {
      u_use = u(i);
    }
    // Calculating actual hazard rate and probability.
    if (re_scale_id == 0){
      haz = exp(log(base_haz + 1e-12) + u_use) + DBL_MIN;
      prob = haz_to_prob(haz);
    } else if (re_scale_id == 1){
      prob = invlogit(logit(base_prob + 1e-12) + u_use) + DBL_MIN;
      haz = prob_to_haz(prob);
    }
    // Running calculation of overall probability of nondetection.
    p_total_evade *= 1 - prob;
  }
  // If we're dealing with a binomial response, then this is only per session.
  if (resp_id == 0){
    p_total_evade = pow(p_total_evade, resp_pars(0));
  }
  if (cov_id == 3){
    f -= dnorm(u(0), Type(0), cov_pars(0), true);
  } else {
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
	    // Total dependence (individual-level random effect).
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
  // Contribution from probability of detection.
  f -= log(1 - p_total_evade);
  return f;
}
