// Template to calculate detection probability at a particular mask
// point for a model with independent random effects.
#include <TMB.hpp>
#include <fenv.h>
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Distance from mask point to all traps.
  DATA_VECTOR(dists);
  // Total number of traps.
  DATA_INTEGER(n_traps);
  // Indicator for dependence structure.
  DATA_INTEGER(dep_id);
  // Parameter values.
  DATA_SCALAR(lambda0);
  DATA_SCALAR(sigma);
  // Parameters of latent variable process.
  DATA_SCALAR(sigma_u);
  // Latent variables.
  PARAMETER_VECTOR(u);
  // Overall probability of nondetection.
  Type p_total_evade = 1;
  // Negative-log joint density of probability of detection and latent variables.
  Type f = 0;
  // Probability of capture.
  for (int i = 0; i < n_traps; i++){
    // Calculating encounter rate.
    Type er = exp(log(lambda0*exp(-pow(dists(i), 2)/(2*pow(sigma, 2)))) + u(i));
    // Calculating probability of detection.
    Type p_detected = 1 - exp(-er);
    // Running calculation of overall probability of nondetection.
    p_total_evade *= 1 - p_detected;
   
    //f -= dnorm(u(i), Type(0), sigma_u, true);
  }
  // Variance-covariance matrix for latent variables.
  matrix<Type> sigma_u_mat(n_traps, n_traps);
  for (int j = 0; j < n_traps; j++){
    for (int k = j; k < n_traps; k++){ 
      if (j == k){
	sigma_u_mat(j, k) = pow(sigma_u, 2);
      } else {
	// Covariance function in here.
	if (dep_id == 0){
	  sigma_u_mat(j, k) = 0;
	  sigma_u_mat(k, j) = 0;
	} else if (dep_id == 1){
	  // Exponential covariance function in here.
	} else if (dep_id == 2){
	  // Matern covariance function in here.
	} else if (dep_id == 3){
	  sigma_u_mat(j, k) = pow(sigma_u, 2);
	  sigma_u_mat(k, j) = pow(sigma_u, 2);
	}
      }
    }
  }
  // Contribution from latent variables (note MVNORM returns the
  // negative-log of the density).
  f += MVNORM(sigma_u_mat)(u);
  // Contribution from probability of detection.
  f -= log(1 - p_total_evade);
  return f;
}
