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
  // Indicator for response type.
  DATA_INTEGER(resp_id);
  // Additional response parameters.
  DATA_VECTOR(resp_pars)
  // Indicator for dependence structure.
  DATA_INTEGER(cov_id);
  // Detection function parmaters.
  DATA_VECTOR(det_pars);
  // Covariance parameters.
  DATA_VECTOR(cov_pars);
  // Latent variables.
  PARAMETER_VECTOR(u);
  // Overall probability of nondetection.
  Type p_total_evade = 1;
  // Negative-log joint density of probability of detection and latent variables.
  Type f = 0;
  // Probability of capture.
  for (int i = 0; i < n_traps; i++){
    // Calculating encounter rate.
    Type er = exp(log(detfn(mask_dists(i), det_pars, detfn_id)) + u(i));
    // Calculating probability of detection.
    Type p_detected = 1 - exp(-er);
    // Running calculation of overall probability of nondetection.
    p_total_evade *= 1 - p_detected;
  }
  // If we're dealing with a binomial response, then this is only per session.
  if (resp_id == 0){
    p_total_evade = pow(p_total_evade, resp_pars(0));
  }
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
	}
      }
    }
  }
  // for (int i = 0; i < 3; i++){
  //   for (int j = 0; j < 3; j++){
  //     std::cout << sigma_u_mat(i, j) << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // exit(1234);
  // Contribution from latent variables (note MVNORM returns the
  // negative-log of the density).
  f += MVNORM(sigma_u_mat)(u);
  // Contribution from probability of detection.
  f -= log(1 - p_total_evade);
  return f;
}
