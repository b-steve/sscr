#ifndef UTILITIES_H
#define UTILITIES_H

template<class Type>
Type dbinom_sscr (const Type &k, const Type &size, const Type &prob, const int &give_log){
  Type out;
  out = exp(lgamma(size + 1) - lgamma(k + 1) - lgamma(size - k + 1))*pow(prob, k)*pow(1 - prob, size - k);
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
  int nr = haz.rows();
  int nc = haz.cols();
  matrix<Type> prob(nr, nc);
  for (int i = 0; i < nr; i++){
    for (int j = 0; j < nc; j++){
      prob(i, j) = haz_to_prob(haz(i, j));
    }
  }
}

// CONVERSION OF PROBABILITIES TO HAZARDS IN HERE.

#endif
