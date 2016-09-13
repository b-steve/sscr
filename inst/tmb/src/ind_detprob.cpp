// Template to calculate detection probability at a particular mask
// point for a model with independent random effects.
#include <TMB.hpp>
#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Distance from mask point to all traps.
  DATA_VECTOR(dists);
  // Total number of traps.
  DATA_INTEGER(n_traps);
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
    // Contribution to joint density from latent variables.
    f -= dnorm(u(i), Type(0), sigma_u, true);
  }
  // Contribution from probability of detection.
  f -= log(1 - p_total_evade);
  return f;
}
