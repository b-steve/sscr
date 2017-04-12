#ifndef UTILITIES_H
#define UTILITIES_H

// Stable binomial PMF.
template<class Type>
Type dbinom_sscr (const Type &k, const Type &size, const Type &prob, const int &give_log){
  Type out;
  out = exp(lgamma(size + 1) - lgamma(k + 1) - lgamma(size - k + 1))*pow(prob, k)*pow(1 - prob, size - k);
  if (give_log){
    out = log(out);
  }
  return out;
}

// Stable Poisson PMF.
template<class Type>
Type dpois_sscr (const Type &x, const Type &lambda, const int &give_log){
  Type out;
  out = pow(lambda, x)*exp(-lambda)/exp(lgamma(x + 1));
  if (give_log){
    out = log(out);
  }
  return out;
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
  return -log(1 - prob);
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

#endif
