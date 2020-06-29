namespace kutils {
  
#define _eps 1e-5 // An alternative limit argument for the first-order IGRMF

// logit
template <class Type>
Type logit(Type p) {
  if ((p < 0) | (p > 1)) {
    Rcout << "p must be in [0, 1]";
    return 0;
  }
  return log(p/(1-p));
}

// inverse logit
template <class Type>
Type logit_inv(Type x) {
  return 1/(1+exp(-x));
}


// soft constraint to zero
template <class Type>
Type soft_zero_sum(vector<Type> x) {
    return dnorm(sum(x), Type(0.0), Type(0.001) * x.size(), true);
}

// Log-logistic distribution
// 
// survival function
template <class Type>
Type St_llogis(Type t, Type alpha, Type lambda) {
  return Type(1.0) / (Type(1.0) + pow(lambda * t, alpha));
}

// density function
template <class Type>
Type ft_llogis(Type t, Type alpha, Type lambda) {
  Type num = lambda * alpha * pow(lambda*t, alpha-Type(1.0));
  Type St  = St_llogis(t, alpha, lambda);
  return num * pow(St, 2);
}

// Mean function
template <class Type>
Type mu_llogis(Type alpha, Type lambda) {
  Type t1 = Type(1.0)/lambda;
  Type t2 = (M_PI/alpha)/sin(M_PI/alpha);
  return t1*t2;
}

// penalized-log precision density
template <class Type>
Type pc_prec(Type x, Type u = Type(1.0), Type alpha = Type(0.01)) {
  Type lambda = -log(alpha) / u;
  Type out = log(0.5 * lambda) - lambda * exp(-0.5 * x) - 0.5 * x;
  return out;
}

// IID
template <class Type>
Type iid(vector<Type> x, Type tau, bool sum_c=true) {
  Type o = Type(0.0);
  o -= pc_prec(tau, Type(1));
  if (sum_c) 
    o -= soft_zero_sum(x);
  o -= dnorm(x, Type(0.0), pow(tau, -0.5), true).sum();
  return o;
}

// prepare Q
template <class Type>
Eigen::SparseMatrix<Type> prepare_Q(matrix<Type> R, Type tau) {
  R = tau * R.array();
  R.diagonal().array()+=_eps;
  return tmbutils::asSparseMatrix(R);
}

// RW1
template <class Type>
Type rw1(vector<Type> x, vector<Type> xid, matrix<Type> R, Type tau,
           bool slope_c=false, bool sum_c=true) {
  Type rwll = Type(0.0);
  if (sum_c) 
    rwll -= soft_zero_sum(x);
  if (slope_c) 
    rwll -= soft_zero_sum(xid*x);
  rwll += density::GMRF(prepare_Q(R, tau))(x);
  return rwll;
}

// RW2
template <class Type>
Type rw2(vector<Type> x, vector<Type> xid, matrix<Type> R, Type tau,
           bool slope_c=true, bool sum_c=true) {
  Type rwll = Type(0.0);
  if (sum_c) 
    rwll -= soft_zero_sum(x);
  if (slope_c) 
    rwll -= soft_zero_sum(xid*x);
  rwll += density::GMRF(prepare_Q(R, tau))(x);
  return rwll;
}

// BESAG
template <class Type>
Type besag(matrix<Type> R, Type tau, vector<Type> u_vec) {
  Type dll = 0.0;
  dll -= pc_prec(tau); // TODO: allow to specify this
  dll += density::GMRF(prepare_Q(R, tau))(u_vec);
  dll -= soft_zero_sum(u_vec);
  return dll;
}

// BYM Model
template <class Type>
Type bym(matrix<Type> R, Type tau_u, Type tau_v, 
  vector<Type> u_vec, vector<Type> v_vec) {
  Type dll = 0.0;
  dll += iid(v_vec, tau_v);
  dll += besag(R, tau_u, u_vec);
  return dll;
}

// DEAN Model: TODO find an example and compare
template <class Type>
Type dean(matrix<Type> R, Type tau, Type phi, 
  vector<Type> u_vec, vector<Type> v_vec) {
  Type dll = 0.0;
  dll += iid(v_vec, tau/(1-phi));
  dll += besag(R, tau/phi, u_vec);
  return dll;
}

// Leroux model: TODO find an example and compare
template <class Type>
Type leroux(matrix<Type> R, Type tau, Type phi, 
  vector<Type> u_vec, vector<Type> v_vec) {
  Type dll = 0.0;
  dll += iid(v_vec, tau*(1-phi));
  dll += besag(R, tau*phi);
  return dll;
}

// BYM2

} // End namespace