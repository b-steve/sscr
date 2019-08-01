// Template to calculate detection probability at a particular mask
// point for a model with some covariance function for the random
// effects.
#include <TMB.hpp>
#include <fenv.h>
#include "utilities.h"
#include "detfns.h"
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
  // Indicator for response type.
  DATA_INTEGER(resp_id);
  // Indicator for dependence structure.
  DATA_INTEGER(cov_id);
  // Indicator for random effect multiplier.
  DATA_INTEGER(mult_id);
  // Indicators for detection parameter link functions.
  DATA_IVECTOR(link_det_ids);
  // Number of detection function parameters.
  int n_det_pars = link_det_ids.size();
  // Indicators for covariance parameter link functions.
  DATA_IVECTOR(link_cov_ids);
  // Number of covariance function parameters.
  int n_cov_pars = link_cov_ids.size();
  // Indicators for additional response parameter link functions.
  DATA_IVECTOR(link_resp_ids);
  // Number of additional response parameters.
  int n_resp_pars = link_resp_ids.size();
  // Numerical offset for the final log().
  DATA_SCALAR(log_offset);
  // Detection function parmaters.
  PARAMETER_VECTOR(link_det_pars);
  // Covariance parameters.
  PARAMETER_VECTOR(link_cov_pars);
  // Additional response parameters.
  PARAMETER_VECTOR(link_resp_pars)
  // Time-of-arrival parameter.
  PARAMETER(link_sigma_toa);
  // Density parameter.
  PARAMETER(link_D);
  // Latent variables.
  PARAMETER_VECTOR(u);
  // Setting a minimum value.
  double dbl_min = 1e-50;
  // Back-transforming detection function parameters.
  vector<Type> det_pars(n_det_pars);
  for (int i = 0; i < n_det_pars; i++){
    if (link_det_ids(i) == 0){
      det_pars(i) = exp(link_det_pars(i));
    } else if (link_det_ids(i) == 1){
      det_pars(i) = 1/(1 + exp(-link_det_pars(i)));
    } else if (link_det_ids(i) == 2){
      det_pars(i) = link_det_pars(i);
    }
  }
  // Back-transforming covariance function parameters.
  vector<Type> cov_pars(n_cov_pars);
  for (int i = 0; i < n_cov_pars; i++){
    if (link_cov_ids(i) == 0){
      cov_pars(i) = exp(link_cov_pars(i));
    } else if (link_cov_ids(i) == 1){
      cov_pars(i) = 1/(1 + exp(-link_cov_pars(i)));
    } else if (link_cov_ids(i) == 2){
      cov_pars(i) = link_cov_pars(i);
    }
  }
  // Back-transforming additional response parameters.
  vector<Type> resp_pars(n_resp_pars);
  for (int i = 0; i < n_resp_pars; i++){
    if (link_resp_ids(i) == 0){
      resp_pars(i) = exp(link_resp_pars(i));
    } else if (link_resp_ids(i) == 1){
      resp_pars(i) = 1/(1 + exp(-link_resp_pars(i)));
    } else if (link_resp_ids(i) == 2){
      resp_pars(i) = link_resp_pars(i);
    }
  }
  // Back-transforming sigma_toa, including readjustment to ms.
  Type sigma_toa = exp(link_sigma_toa)/1000;
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
    // Calculating baseline hazard rate and probability. Needs to be
    // done separately for CMP distribution.
    if (resp_id == 2){
      base_haz = prob_to_haz(detfn(mask_dists(i), det_pars, detfn_id));
      base_prob = cmp_haz_to_prob(base_haz, resp_pars(0));
    } else {
      base_prob = detfn(mask_dists(i), det_pars, detfn_id);
      base_haz = prob_to_haz(base_prob);
    }
    if (cov_id == 6){
      u_use = 0;
    } else if (cov_id == 3){
      u_use = u(0);
    } else {
      u_use = u(i);
    }
    // Calculating actual hazard rate and probability.
    if (mult_id == 0){
      haz = base_haz*exp(u_use) + dbl_min;
    } else if (mult_id == 1){
      haz = base_prob*exp(u_use) + dbl_min;
    }
    // Converting to a probability.
    if (resp_id == 2){
      prob = cmp_haz_to_prob(haz, resp_pars(0));
    } else {
      prob = haz_to_prob(haz);
    }
    // Running calculation of overall probability of nondetection.
    p_total_evade *= 1 - prob;
  }
  // If we're dealing with a binomial response, then this is only per session.
  if (resp_id == 0){
    p_total_evade = pow(p_total_evade, resp_pars(0));
  }
  if (cov_id == 3){
    f -= dnorm(u(0), cov_pars(0), cov_pars(1), true);
  } else if (cov_id != 6){
    // Variance-covariance matrix for latent variables.
    matrix<Type> sigma_u_mat(n_traps, n_traps);
    for (int j = 0; j < n_traps; j++){
      for (int k = j; k < n_traps; k++){ 
	if (j == k){
	sigma_u_mat(j, k) = pow(cov_pars(1), 2);
	} else {
	  if (cov_id == 0){
	    // Independent random effects
	    sigma_u_mat(j, k) = 0;
	    sigma_u_mat(k, j) = 0;
	  } else if (cov_id == 1){
	    // Exponential covariance function.
	    sigma_u_mat(j, k) = pow(cov_pars(1), 2)*exp(-trap_dists(j, k)/cov_pars(2));
	    sigma_u_mat(k, j) = pow(cov_pars(1), 2)*exp(-trap_dists(j, k)/cov_pars(2));
	  } else if (cov_id == 2){
	    // Matern covariance function.
	  } else if (cov_id == 3){
	    // Total dependence (individual-level random effect).
	    sigma_u_mat(j, k) = pow(cov_pars(1), 2);
	    sigma_u_mat(k, j) = pow(cov_pars(1), 2);
	  } else if (cov_id == 4){
	    // Linear combination of exponential covariance functions.
	    sigma_u_mat(j, k) = pow(cov_pars(1), 2)*(cov_pars(2)*exp(-cov_pars(3)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(3)*trap_dists(j, k)), 2) + (1 - cov_pars(2))*exp(-cov_pars(4)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(4)*trap_dists(j, k)), 2));
	    sigma_u_mat(k, j) = pow(cov_pars(1), 2)*(cov_pars(2)*exp(-cov_pars(3)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(3)*trap_dists(j, k)), 2) + (1 - cov_pars(2))*exp(-cov_pars(4)*trap_dists(j, k))*1/pow(1 + exp(-cov_pars(4)*trap_dists(j, k)), 2));
	  } else if (cov_id == 5){
	    // Squared exponential covariance function.
	    sigma_u_mat(j, k) = pow(cov_pars(1), 2)*exp(-pow(trap_dists(j, k), 2)/pow(cov_pars(2), 2));
	    sigma_u_mat(k, j) = pow(cov_pars(1), 2)*exp(-pow(trap_dists(j, k), 2)/pow(cov_pars(2), 2));
	  } else {
	    exit(1111);
	  }
	}
      }
    }
    // Contribution from latent variables (note MVNORM returns the
    // negative-log of the density). Subtracting off mu.u because
    // MVNORM assumes zero mean.
    f += MVNORM(sigma_u_mat)(u - cov_pars(0));
  }
  // Contribution from probability of detection.
  f -= log(1 - p_total_evade + log_offset);
  return f;
}
