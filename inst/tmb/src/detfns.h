#ifndef DETFNS_H
#define DETFNS_H

template<class Type>
Type detfn_hn (const Type &d, const vector<Type> &det_pars)
{
  return det_pars(0)*exp(-pow(d, 2)/(2*pow(det_pars(1), 2)));
}

template<class Type>
Type detfn_hr (const Type &d, const vector<Type> &det_pars)
{
  return det_pars(0)*(1 - exp(-pow(d/det_pars(1), -det_pars(2))));
}

template<class Type>
Type detfn_hhn (const Type &d, const vector<Type> &det_pars)
{
  return haz_to_prob(detfn_hn(d, det_pars));
}

template<class Type>
Type detfn_hhr (const Type &d, const vector<Type> &det_pars)
{
  return haz_to_prob(det_pars(0)*(1 - exp(-pow(d/det_pars(1), -det_pars(2)))));
}

template<class Type>
Type detfn (const Type &d, const vector<Type> &det_pars, const int &detfn_id)
{
  Type out;
  if (detfn_id == 0){
    out = detfn_hn(d, det_pars);
  } else if (detfn_id == 1){
    out = detfn_hr(d, det_pars);
  } else if (detfn_id == 2){
    out = detfn_hhn(d, det_pars);
  } else if (detfn_id == 3){
    out = detfn_hhr(d, det_pars);
  }
  return out;
}

#endif
