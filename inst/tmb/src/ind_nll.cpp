// Template to calculate the negative log-likelihood for a model with
// independent random effects.
#include <TMB.hpp>
#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Reading in data.
  DATA_MATRIX(capt);
  DATA_MATRIX(mask_dists);
  DATA_INTEGER(n);
  DATA_INTEGER(n_traps);
  DATA_INTEGER(n_mask);
  DATA_SCALAR(mask_area);
  // Detection probabilities (from laplace_probs.cpp).
  DATA_VECTOR(prob_det);
  // Declaring parameters.
  DATA_SCALAR(log_lambda0);
  DATA_SCALAR(log_sigma);
  DATA_SCALAR(log_sigma_u);
  PARAMETER_MATRIX(u);
  Type lambda0 = exp(log_lambda0);
  ADREPORT(lambda0);
  Type sigma = exp(log_sigma);
  ADREPORT(sigma);
  Type sigma_u = exp(log_sigma_u);
  ADREPORT(sigma_u);
  // Hazard rates for mask/trap combinations.
  matrix<Type> haz_mat(n_mask, n_traps);
  // The sum of mask probabilities.
  Type sum_prob_det = 0;
  for (int i = 0; i < n_mask; i++){
    Type p_undet = Type(1);
    for (int j = 0; j < n_traps; j++){
      haz_mat(i, j) = lambda0*exp(-pow(mask_dists(i, j), 2)/(2*pow(sigma, 2)));
    }
    sum_prob_det += prob_det(i);
  }
  // PMF for activity centres across the mask.
  vector<Type> f_loc(n_mask);
  f_loc = prob_det/sum_prob_det;
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
	//Type e_count = exp(log(haz_mat(j, k)));
	integrand_mask *= dpois(capt(i, k), e_count, false);
      }
      integrand += integrand_mask;
    }
    //integrand *= mask_area;
    log_sum_integrands += log(integrand + DBL_MIN);
  }
  f -= log_sum_integrands;
  Type esa = mask_area*sum_prob_det;
  ADREPORT(esa);
  Type D = n/esa;
  ADREPORT(D);
  // Extra bit that falls out of log-likelihood.
  f -= -n*log(sum_prob_det);
  // Comment out for conditional likelihood.
  //f -= dpois(Type(n), D*esa, true);
  // Contribution from latent variables.
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n_traps; j++){
      f -= dnorm(u(i, j), Type(0), sigma_u, true);
    }
  }
  return f;
}