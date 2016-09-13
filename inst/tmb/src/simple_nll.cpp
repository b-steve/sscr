// Template to calculate the negative log-likelihood for a model with
// no random effects.
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
  // Indicator for poisson or binary.
  DATA_INTEGER(resp_id);
  // Declaring parameters.
  PARAMETER(log_lambda0);
  PARAMETER(log_sigma);
  Type lambda0 = exp(log_lambda0);
  ADREPORT(lambda0);
  Type sigma = exp(log_sigma);
  ADREPORT(sigma);
  // Hazard rates for mask/trap combinations.
  matrix<Type> haz_mat(n_mask, n_traps);
  // Detection probabilities for mask/trap combinations.
  matrix<Type> prob_mat(n_mask, n_traps);
  // Detection probabilities for each mask point.
  vector<Type> prob_det(n_mask);
  // The sum of mask probabilities.
  Type sum_prob_det = 0;
  for (int i = 0; i < n_mask; i++){
    Type p_undet = Type(1);
    for (int j = 0; j < n_traps; j++){
      haz_mat(i, j) = lambda0*exp(-pow(mask_dists(i, j), 2)/(2*pow(sigma, 2)));
      prob_mat(i, j) = 1 - exp(-haz_mat(i, j));
      p_undet *= 1 - prob_mat(i, j);
    }
    prob_det(i) = 1 - p_undet;
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
      Type integrand_mask = 1;
      for (int k = 0; k < n_traps; k++){
	if (resp_id == 0){
	  integrand_mask *= pow(prob_mat(j, k), capt(i, k))*pow(1 - prob_mat(j, k), 1 - capt(i, k));
	  // Function dbinom() doesn't seem to work properly.
	  //integrand_mask *= dbinom(capt(i, k), Type(1), prob_mat(j, k), false);
	} else if (resp_id == 1){
	  integrand_mask *= dpois(capt(i, k), haz_mat(j, k), false);
	}
      }
      integrand += integrand_mask;
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
  // Comment out for conditional likelihood.
  //f -= dpois(Type(n), D*esa, true);
  //std::cout << "D: " << D << ", lambda0: " << lambda0 << ", sigma: " << sigma << ", f: " << f << std::endl;
  return f;
}
