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
  // Indicator for detection function.
  DATA_INTEGER(detfn_id);
  // Indicator for dependence structure.
  DATA_INTEGER(cov_id);
  // Indicator for random effect multiplier.
  DATA_INTEGER(mult_id);
  // Detection probabilities (from cov_detprob.cpp).
  DATA_VECTOR(det_probs);
  // Indicator for time-of-arrival data.
  DATA_INTEGER(toa_id);
  // Time-of-arrival data.
  DATA_MATRIX(toa_ssq);
  // Indicator for conditional likelihood.
  DATA_INTEGER(conditional_n);
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
  PARAMETER_MATRIX(u);
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
    } else if (link_resp_ids(i) == 3){
      resp_pars(i) = exp(link_resp_pars(i)) + 1;
    }
  }
  // Back-transforming sigma_toa, including readjustment to ms.
  Type sigma_toa = exp(link_sigma_toa)/1000;
  // Back-transforming density parameter.
  Type D = exp(link_D);
  // Hazard rates for mask/trap combinations.
  matrix<Type> haz_mat(n_mask, n_traps);
  // Detection probabilities for mask/trap combinations.
  matrix<Type> prob_mat(n_mask, n_traps);
  // Needs to be done separately for non-binomial/non-Poisson distributions.
  if (resp_id == 2){
    for (int i = 0; i < n_mask; i++){
      for (int j = 0; j < n_traps; j++){
	haz_mat(i, j) = prob_to_haz(detfn(mask_dists(i, j), det_pars, detfn_id));
      }
    }
    prob_mat = cmp_haz_to_prob(haz_mat, resp_pars(0));
  } else if (resp_id == 3){
    for (int i = 0; i < n_mask; i++){
      for (int j = 0; j < n_traps; j++){
	haz_mat(i, j) = prob_to_haz(detfn(mask_dists(i, j), det_pars, detfn_id));
      }
    }
    prob_mat = nb_haz_to_prob(haz_mat, resp_pars(0));
  } else if (resp_id == 4){
    for (int i = 0; i < n_mask; i++){
      for (int j = 0; j < n_traps; j++){
	haz_mat(i, j) = prob_to_haz(detfn(mask_dists(i, j), det_pars, detfn_id));
      }
    }
    prob_mat = nba_haz_to_prob(haz_mat, resp_pars(0));
  } else {
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
  // Declaring the vector of random effects to use.
  Type u_use;
  // Declaring expected encounter-rate and detection probabilities.
  Type e_count;
  Type e_prob;
  // Joint density of data and latent variables.
  Type f = 0;
  // Likelihood contributions from capture histories.
  Type log_sum_integrands = 0;
  for (int i = 0; i < n; i++){
    Type integrand = 0;
    for (int j = 0; j < n_mask; j++){
      Type integrand_mask = 0;
      for (int k = 0; k < n_traps; k++){
	if (cov_id == 6){
	  u_use = 0;
	} else if (cov_id == 3){
	  u_use = u(i, 0);
	} else {
	  u_use = u(i, k);
	}
	// Expected counts and probabilities.
	if (mult_id == 0){
	  e_count = haz_mat(j, k)*exp(u_use) + dbl_min;
	} else if (mult_id == 1){
	  e_count = prob_mat(j, k)*exp(u_use) + dbl_min;
	}
	if (resp_id == 0){
	  e_prob = haz_to_prob(e_count);
	  integrand_mask += dbinom_sscr(capt(i, k), resp_pars(0), e_prob, true);
	} else if (resp_id == 1){
	  integrand_mask += dpois_sscr(capt(i, k), e_count, true);
	} else if (resp_id == 2){
	  integrand_mask += dcompois2(capt(i, k), e_count, resp_pars(0), true);
	} else if (resp_id == 3){
	  integrand_mask += dnbinom_sscr(capt(i, k), e_count, resp_pars(0), true);
	} else if (resp_id == 4){
	  integrand_mask += dnbinomalpha_sscr(capt(i, k), e_count, resp_pars(0), true);
	}
      }
      // Time-of-arrival component.
      if (toa_id == 1){
	integrand_mask += (1 - n_dets(i))*log(sigma_toa) - (toa_ssq(i, j)/(2*pow(sigma_toa, 2)));
      }
      integrand += exp(integrand_mask);
    }
    log_sum_integrands += log(integrand + dbl_min);
  }
  f -= log_sum_integrands;
  // Extra bit that falls out of log-likelihood.
  f -= -n*log(sum_det_probs);
  // Likelihood component due to n.
  if (conditional_n == 0){
    f -= dpois(Type(n), Type(D*mask_area*sum_det_probs), true);
  }
  if (cov_id == 3){
    for (int i = 0; i < n; i++){
      f -= dnorm(u(i, 0), cov_pars(0), cov_pars(1), true);
    }
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
	    // Total dependence (individual-level random effects).
	    // NOTE: This is actually handled above.
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
    for (int i = 0; i < n; i++){
      // Contribution from latent variables (note MVNORM returns the
      // negative-log of the density).
      vector<Type> v(u.cols());
      for (int ii = 0; ii < u.cols(); ii++){
	v(ii) = u.row(i)(ii) - cov_pars(0);
      }
      f += MVNORM(sigma_u_mat)(v);
    }
  }
  return f;
}
