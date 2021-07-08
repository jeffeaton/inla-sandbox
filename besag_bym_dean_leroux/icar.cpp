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
  PARAMETER(log_icar_tau);
  Type icar_tau = exp(log_icar_tau);
  dll += kutils::besag(R_, icar_tau, icar_vec);

  PARAMETER(intercept);
  dll -= dnorm(intercept, Type(0), Type(10), true);

  PARAMETER(log_sd_gauss);
  Type sd_gauss = exp(log_sd_gauss);
  dll -= kutils::pc_prec(pow(sd_gauss, -2));

  // Data likelihood
  for (int i=0; i < rent.size(); i++) {
    Type eta = intercept + icar_vec(location(i));
    dll -= dnorm(rent(i), eta, sd_gauss, true);
  }

  Type t_gau = 1/pow(sd_gauss,2);
  Type t_loc = icar_tau;

  ADREPORT(t_gau);
  ADREPORT(t_loc);
  return dll;
}
