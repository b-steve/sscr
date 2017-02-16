#ifndef UTILITIES_H
#define UTILITIES_H

#endif

template<class Type>
Type dbinom_sscr (const Type &k, const Type &size, const Type &prob, const int &give_log){
  Type out;
  out = exp(lgamma(size + 1) - lgamma(k + 1) - lgamma(size - k + 1))*pow(prob, k)*pow(1 - prob, size - k);
  if (give_log){
    out = log(out);
  }
  return out;
}
