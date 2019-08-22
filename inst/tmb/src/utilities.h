#ifndef UTILITIES_H
#define UTILITIES_H

// Stable binomial PMF.
template<class Type>
Type dbinom_sscr (const Type &k, const Type &size, const Type &prob, const int &give_log){
  Type out;
  out = exp(lgamma(size + 1) - lgamma(k + 1) - lgamma(size - k + 1))*pow(prob, k)*pow(1 - prob, size - k);
  if (give_log){
    out = log(out + DBL_MIN);
  }
  return out;
}

// Stable Poisson PMF.
template<class Type>
Type dpois_sscr (const Type &x, const Type &lambda, const int &give_log){
  Type out;
  out = pow(lambda, x)*exp(-lambda)/exp(lgamma(x + 1));
  if (give_log){
    out = log(out + DBL_MIN);
  }
  return out;
}
template<class Type>
Type dpois_sscr (const int &x, const Type &lambda, const int &give_log){
  double d_x = x;
  Type out;
  out = pow(lambda, d_x)*exp(-lambda)/exp(lgamma(d_x + 1));
  if (give_log){
    out = log(out + DBL_MIN);
  }
  return out;
}

// Negative binomial distribution.
template<class Type>
Type dnbinom_sscr(const int &x, const Type &mu, const Type &size, const int &give_log){
  Type x_var = mu + pow(mu, 2)/size;
  dbinom_robust(x, log(mu + DBL_MIN), log(x_var - mu + DBL_MIN), give_log);
}

// Converting hazards to probabilities.

// For scalars.
template<class Type>
Type haz_to_prob (const Type &haz){
  return 1 - exp(-haz);
}

// For vectors.
template<class Type>
vector<Type> haz_to_prob (const vector<Type> &haz){
  int n = haz.size();
  vector<Type> prob(n);
  for (int i = 0; i < n; i++){
    prob(i) = haz_to_prob(haz(i));
  }
  return prob;
}

// For matrices.
template<class Type>
matrix<Type> haz_to_prob (const matrix<Type> &haz){
  int nr = haz.col(1).size();
  int nc = haz.row(1).size();
  matrix<Type> prob(nr, nc);
  for (int i = 0; i < nr; i++){
    for (int j = 0; j < nc; j++){
      prob(i, j) = haz_to_prob(haz(i, j));
    }
  }
  return prob;
}

// Converting probabilities to hazards.
// For scalars.
template<class Type>
Type prob_to_haz (const Type &prob){
  return -log(1 - prob + DBL_MIN);
}

// For vectors.
template<class Type>
vector<Type> prob_to_haz (const vector<Type> &prob){
  int n = prob.size();
  vector<Type> haz(n);
  for (int i = 0; i < n; i++){
    haz(i) = prob_to_haz(prob(i));
  }
  return haz;
}

// For matrices.
template<class Type>
matrix<Type> prob_to_haz (const matrix<Type> &prob){
  int nr = prob.col(1).size();
  int nc = prob.row(1).size();
  matrix<Type> haz(nr, nc);
  for (int i = 0; i < nr; i++){
    for (int j = 0; j < nc; j++){
      haz(i, j) = prob_to_haz(prob(i, j));
    }
  }
  return haz;
}

// Calculating detection probability for CMP distribution.
// For scalars.
template<class Type>
Type cmp_haz_to_prob(const Type &haz, const Type &nu){
  Type loglambda = compois_calc_loglambda(log(haz + DBL_MIN), nu);
  return 1 - 1/exp(compois_calc_logZ(loglambda, nu));
}

// For vectors.
template<class Type>
vector<Type> cmp_haz_to_prob (const vector<Type> &haz, const Type &nu){
  int n = haz.size();
  vector<Type> prob(n);
  for (int i = 0; i < n; i++){
    prob(i) = cmp_haz_to_prob(haz(i), nu);
  }
  return prob;
}

// For matrices.
template<class Type>
matrix<Type> cmp_haz_to_prob (const matrix<Type> &haz, const Type &nu){
  int nr = haz.col(1).size();
  int nc = haz.row(1).size();
  matrix<Type> prob(nr, nc);
  for (int i = 0; i < nr; i++){
    for (int j = 0; j < nc; j++){
      prob(i, j) = cmp_haz_to_prob(haz(i, j), nu);
    }
  }
  return prob;
}


// Calculating detection probability for NB distribution.
// For scalars.
template<class Type>
Type nb_haz_to_prob(const Type &haz, const Type &size){
  return 1 - pow((size/(size + haz)), size);
}

// For vectors.
template<class Type>
vector<Type> nb_haz_to_prob (const vector<Type> &haz, const Type &size){
  int n = haz.size();
  vector<Type> prob(n);
  for (int i = 0; i < n; i++){
    prob(i) = nb_haz_to_prob(haz(i), size);
  }
  return prob;
}

// For matrices.
template<class Type>
matrix<Type> nb_haz_to_prob (const matrix<Type> &haz, const Type &size){
  int nr = haz.col(1).size();
  int nc = haz.row(1).size();
  matrix<Type> prob(nr, nc);
  for (int i = 0; i < nr; i++){
    for (int j = 0; j < nc; j++){
      prob(i, j) = nb_haz_to_prob(haz(i, j), size);
    }
  }
  return prob;
}


/* // Logistic function. */

/* // For scalars. */

/* template<class Type> */
/* Type logit (const Type &prob){ */
/*   return log(prob/(1 - prob)); */
/* } */

/* // For vectors. */
/* template<class Type> */
/* vector<Type> logit (const vector<Type> &prob){ */
/*   int n = prob.size(); */
/*   vector<Type> logodds(n); */
/*   for (int i = 0; i < n; i++){ */
/*     logodds(i) = logit(prob(i)); */
/*   } */
/*   return logodds; */
/* } */

/* // For matrices. */
/* template<class Type> */
/* matrix<Type> logit (const matrix<Type> &prob){ */
/*   int nr = prob.col(1).size(); */
/*   int nc = prob.row(1).size(); */
/*   matrix<Type> logodds(nr, nc); */
/*   for (int i = 0; i < nr; i++){ */
/*     for (int j = 0; j < nc; j++){ */
/*       logodds(i, j) = logit(prob(i, j)); */
/*     } */
/*   } */
/*   return logodds; */
/* } */

/* // Inverse logistic function. */

/* template<class Type> */
/* Type invlogit (const Type &logodds){ */
/*   return 1/(1 + exp(-logodds)); */
/* } */

/* // For vectors. */
/* template<class Type> */
/* vector<Type> invlogit (const vector<Type> &logodds){ */
/*   int n = logodds.size(); */
/*   vector<Type> prob(n); */
/*   for (int i = 0; i < n; i++){ */
/*     prob(i) = invlogit(logodds(i)); */
/*   } */
/*   return prob; */
/* } */

/* // For matrices. */
/* template<class Type> */
/* matrix<Type> invlogit (const matrix<Type> &logodds){ */
/*   int nr = logodds.col(1).size(); */
/*   int nc = logodds.row(1).size(); */
/*   matrix<Type> prob(nr, nc); */
/*   for (int i = 0; i < nr; i++){ */
/*     for (int j = 0; j < nc; j++){ */
/*       prob(i, j) = invlogit(logodds(i, j)); */
/*     } */
/*   } */
/*   return prob; */
/* } */


#endif
