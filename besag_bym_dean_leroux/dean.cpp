
#include <TMB.hpp>
#include <fenv.h> // Extra line needed
#include "k_utils.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  feraiseexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);
  parallel_accumulator<Type> dll(this); // openMP only works on Linux

  // survival data
  DATA_VECTOR(rent);

  // ICAR
  DATA_IVECTOR(location);
  DATA_MATRIX(R_);

  PARAMETER_VECTOR(icar_vec);
  PARAMETER_VECTOR(iid_vec);
  PARAMETER(log_tau);
  PARAMETER(logit_phi);
  Type tau = exp(log_tau), phi = kutils::logit_inv(logit_phi);
  dll += kutils::dean(R_, tau, phi, icar_vec, iid_vec);

  PARAMETER(intercept);
  dll -= dnorm(intercept, Type(0), Type(10), true);

  PARAMETER(log_sd_gauss);
  Type sd_gauss = exp(log_sd_gauss);
  dll -= kutils::pc_prec(pow(sd_gauss, -2));

  // Data likelihood
  for (int i=0; i < rent.size(); i++) {
    Type eta = intercept + icar_vec(location(i)) + iid_vec(location(i));
    dll -= dnorm(rent(i), eta, sd_gauss, true);
  }

  Type t_gau = 1/pow(sd_gauss,2);
  Type t_icar = tau / phi;
  Type t_iid = tau / (1 - phi);
  vector<Type> uv = icar_vec + iid_vec;

  ADREPORT(t_gau);
  ADREPORT(t_icar);
  ADREPORT(t_iid);
  ADREPORT(tau);
  ADREPORT(phi);
  ADREPORT(uv);
  return dll;
}