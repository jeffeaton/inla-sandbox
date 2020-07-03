#' Data and example taken from: https://becarioprecario.bitbucket.io/inla-gitbook/ch-temporal.html#sec:spacetime
#'

#' # Setup
#'
#' ## Load packages
library(INLA)
library(DClusterm)
library(tidyverse)
data(brainNM)

#' ## Utility functions
#'
#' @param code model code as a character string
#'
#' @return name of the loaded DLL
#' 
tmb_compile_and_load <- function(code) {
  f <- tempfile(fileext = ".cpp")
  writeLines(mod, f)
  TMB::compile(f)
  dyn.load(TMB::dynlib(tools::file_path_sans_ext(f)))
  basename(tools::file_path_sans_ext(f))
}


nm.adj <- poly2nb(brainst@sp)
adj.mat <- as(nb2mat(nm.adj, style = "B"), "Matrix")

data <- brainst@data

#' Note that area ID = 11 does not have any observed events and very low expected cases.
#' This creates some issues later for the improper models.
#' 

filter(data, ID == 11)


#' Code factor versions of ID variables

data$ID <- as.integer(data$ID)
data$IDf <- factor(sprintf("%02d", data$ID))
data$Yearf <- factor(data$Year, unique(data$Year))

data$ID.Year <- data$Year - 1973 + 1
data$ID2 <- data$ID

data$area.year <- interaction(data$IDf, data$Year)
data$id.area.year <- as.integer(data$area.year)



# Prior of precision
prec.prior <- list(prec = list(param = c(0.001, 0.001)))
diagval <- INLA:::inla.set.f.default()$diagonal

brain.st <- inla(Observed ~ 1 + f(Year, model = "rw1",
                                  hyper = prec.prior) + 
                   f(as.numeric(ID), model = "besag", graph = adj.mat,
                     hyper = prec.prior),
                 data = data, E = Expected, family = "poisson",
                 control.predictor = list(compute = TRUE, link = 1))

summary(brain.st)

names(inla.models()$group)


#' Use options control.inla = list(strategy = "gaussian", int.strategy = "eb")
#' so that we are doing the same thing as in TMB.
#'
#' INLA defaults
#' - Intercept: flat prior
#' - Fixed effects: N(0, prec = 0.001)

brain.st2 <- inla(Observed ~ 1 + 
                    f(as.numeric(ID2), model = "besag", graph = adj.mat,
                      group = ID.Year, control.group = list(model = "rw1"),
                      hyper = prec.prior),
                  data = data, E = Expected, family = "poisson",
                  control.compute = list(config = TRUE),
                  control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                  control.predictor = list(compute = TRUE, link = 1))


#' ## Model 1: intercept only
#'
#' INLA default is a improper flat prior on the intercept

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  vector<Type> mu(beta0 + log(E));
  val -= dpois(y, exp(mu), true).sum();

  return val;
}
'

dll <- tmb_compile_and_load(mod)

tmbdata <- list(y = data$Observed,
                E = data$Expected)
tmbpar <- list(beta0 = 0)

obj <- TMB::MakeADFun(data = tmbdata,
                     parameters = tmbpar,
                     random = c(),
                     DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)


sdr1 <- TMB::sdreport(obj)

inlafit <- inla(Observed ~ 1,
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"))

mgcvfit1 <- gam(Observed ~ 1, data = data, offset = log(Expected), family = "poisson")

summary(sdr1)

summary(inlafit)
inlafit$summary.fixed





#' ## Model 2: ICAR model

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_space);
  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR area model
  DATA_SCALAR(Qrank);


  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space);
  // Note: dgamma() is parameterised as (shape, scale); INLA parameterised as (shape, rate)
  val -= dlgamma(log_prec_space, Type(0.001), Type(1.0 / 0.001), true);

  PARAMETER_VECTOR(u_space);
  val -=  Qrank * 0.5 * log_prec_space -
             0.5 * exp(log_prec_space) * (u_space * (Q * u_space)).sum();
  val -= dnorm(u_space.sum(), Type(0.0), Type(0.001) * u_space.size(), true); // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_space * u_space +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  return val;
}
'

dll <- tmb_compile_and_load(mod)

Q <- diag(rowSums(adj.mat)) - adj.mat
Qadj <- Q + Matrix::Diagonal(ncol(Q), rep(1e-6, ncol(Q)))

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space = Matrix::sparse.model.matrix(~0 + IDf, data),
                Q = Qadj,
                Qrank = as.integer(rankMatrix(Q)))

tmbpar <- list(beta0 = 0,
               log_prec_space = 0,
               u_space = numeric(ncol(tmbdata$Z_space)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_space"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr,
                 control = list(iter.max = 1000,
                                eval.max = 1000))

tmbfit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

sdr2 <- TMB::sdreport(obj)
summary(sdr2, "all")

prec.prior <- list(prec = list(param = c(0.001, 0.001)))

inlafit <- inla(Observed ~ f(as.integer(IDf), model = "besag",
                             hyper = prec.prior, graph = adj.mat, constr = TRUE,
                             diagonal = diagval),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))


summary(inlafit)
inlafit$internal.summary.hyperpar
inlafit$misc$configs$config[[1]]$theta

summary(sdr, "fixed")

inlafit$summary.fixed
inlafit$summary.random[[1]]

plot(inlafit$summary.random[[1]][ , 2], summary(sdr, "random")[-1, 1])
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][ , 3], summary(sdr, "random")[-1, 2])
abline(0, 1, col = "red")


#' #### INLA with custom Cmatrix

diagval <- INLA:::inla.set.f.default()$diagonal

inlafitC <- inla(Observed ~ f(as.integer(IDf), model = "generic0", Cmatrix = Q,
                              hyper = prec.prior, diagonal = diagval, constr = TRUE),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

#' These match:

summary(inlafit)
summary(inlafitC)

#' By setting `constr = TRUE` above, INLA infers the rank deficiency of one for the
#' `Cmatrix`.

grep("rank", inlafitC$logfile, value = TRUE)

#' #### With `constr = FALSE`
#'
#' If I change `constr = FALSE`, for the ICAR model INLA still calculates the rank
#' deficiency of 1 based on the number of connected components of the graph. But for
#' the `"generic0"` version, INLA no longer detects rank deficiency of the
#' `Cmatrix`.

inlafit_unconstr <- inla(Observed ~ f(as.integer(IDf), model = "besag",
                                      hyper = prec.prior, graph = adj.mat, constr = FALSE,
                                      diagonal = diagval),
                         data = data, E = Expected, family = "poisson",
                         control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                         control.compute = list(config = TRUE))

inlafitC_unconstr_bad <- inla(Observed ~ f(as.integer(IDf), model = "generic0", Cmatrix = Q,
                                           hyper = prec.prior, diagonal = diagval, constr = FALSE),
                              data = data, E = Expected, family = "poisson",
                              control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                              control.compute = list(config = TRUE))

grep("rank", inlafit_unconstr$logfile, value = TRUE)
grep("rank", inlafitC_unconstr_bad$logfile, value = TRUE)

#' Consequently, the precision estimate is different:

summary(inlafit_unconstr)
summary(inlafitC_unconstr_bad)

#' We need to specify `f(..., rankdef = 1)` for the results to align

inlafitC_unconstr <- inla(Observed ~ f(as.integer(IDf), model = "generic0", Cmatrix = Q,
                                       hyper = prec.prior, diagonal = diagval, constr = FALSE, rankdef = 1),
                          data = data, E = Expected, family = "poisson",
                          control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                          control.compute = list(config = TRUE))

grep("rank", inlafitC_unconstr$logfile, value = TRUE)

summary(inlafit_unconstr)
summary(inlafitC_unconstr)

#' #### MRF in `mgcv`

rownames(Q) <- colnames(Q) <- levels(data$IDf)

mgcvfit2 <- gam(Observed ~ 0 + s(IDf, bs = "mrf", xt = list(penalty = as.matrix(Q))),
                data = data, offset = log(Expected), family = "poisson")

sm <- smoothCon(s(IDf, bs = "mrf", xt = list(penalty = as.matrix(Q))), data,
                scale.penalty = FALSE)

sm_abs <- smoothCon(s(IDf, bs = "mrf", xt = list(penalty = as.matrix(Q))), data,
                absorb.cons = TRUE, scale.penalty = FALSE)

vcov(mgcvfit2)

mgcvfit2 <- 


#' ### Model 2b: express as TMB GMRF

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_space);
  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR area model
  DATA_SCALAR(Qrank);


  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space);
  // Note: dgamma() is parameterised as (shape, scale); INLA parameterised as (shape, rate)
  val -= dlgamma(log_prec_space, Type(0.001), Type(1.0 / 0.001), true);

  PARAMETER_VECTOR(u_space);
  val +=  SCALE(GMRF(Q), exp(-0.5 * log_prec_space))(u_space);
  val -= -(Q.cols() - Qrank) * 0.5 * (log_prec_space - log(2 * PI)); // adjust GMRF for rank deficiency
  val -= dnorm(sum(u_space), Type(0.0), Type(0.001) * u_space.size(), true); // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_space * u_space +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  return val;
}
'

dll <- tmb_compile_and_load(mod)

Q <- diag(rowSums(adj.mat)) - adj.mat

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space = Matrix::sparse.model.matrix(~0 + IDf, data),
                Q = Q,
                Qrank = as.integer(rankMatrix(Q)))

tmbpar <- list(beta0 = 0,
               log_prec_space = 0,
               u_space = numeric(ncol(tmbdata$Z_space)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_space"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr,
                 control = list(iter.max = 1000,
                                eval.max = 1000))

tmbfit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

sdr2b <- TMB::sdreport(obj)
summary(sdr2b, "all")


#' Hyper parameters

c(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))
summary(sdr2, "fixed")
summary(sdr2b, "fixed")

#' Intercept

inlafit$summary.fixed[ , 1:2]

summary(sdr2, "random")[1, ]
summary(sdr2b, "random")[1, ]

#' Random effects

inlafit$summary.random[[1]]

plot(summary(sdr2, "random")[-1, 1], summary(sdr2b, "random")[-1, 1])
abline(0, 1, col = "red")

plot(summary(sdr2, "random")[-1, 2], summary(sdr2b, "random")[-1, 2])
abline(0, 1, col = "red")



#' ### Model 2c: non-centered parameterisation

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_space);
  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR area model
  DATA_SCALAR(Qrank);


  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space);
  // Note: TMB dlgamma() is parameterised as (shape, scale); INLA parameterised as (shape, rate)
  val -= dlgamma(log_prec_space, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space(exp(-0.5 * log_prec_space));

  PARAMETER_VECTOR(u_raw_space);
  vector<Type> u_space(u_raw_space * sigma_space);

  val += GMRF(Q)(u_raw_space);
  val -= dnorm(sum(u_raw_space), Type(0.0), Type(0.001) * u_raw_space.size(), true); // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_space * u_space +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

Q <- diag(rowSums(adj.mat)) - adj.mat

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space = Matrix::sparse.model.matrix(~0 + IDf, data),
                Q = Q,
                Qrank = as.integer(rankMatrix(Q)))

tmbpar <- list(beta0 = 0,
               log_prec_space = 0,
               u_raw_space = numeric(ncol(tmbdata$Z_space)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

tmbfit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

sdr2c <- TMB::sdreport(obj)
summary(sdr2c, "all")


#' Hyper parameters

summary(sdr2, "fixed")
summary(sdr2b, "fixed")
summary(sdr2c, "fixed")

#' Intercept

summary(sdr2, "random")[1, ]
summary(sdr2b, "random")[1, ]
summary(sdr2c, "random")[1, ]

#' Random effects

inlafit$summary.random[[1]]

plot(summary(sdr2, "random")[-1, 1], summary(sdr2c, "report")[ , 1])
abline(0, 1, col = "red")

plot(summary(sdr2, "random")[-1 , 2], summary(sdr2c, "report")[ , 2])
abline(0, 1, col = "red")


#' ### Model 2d: explicit sum-to-zero constraint

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(L_space); 
  DATA_SPARSE_MATRIX(Z_space);
  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR area model

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space);
  // Note: TMB dlgamma() is parameterised as (shape, scale); INLA parameterised as (shape, rate)
  val -= dlgamma(log_prec_space, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space(exp(-0.5 * log_prec_space));

  PARAMETER_VECTOR(u_raw_space);
  vector<Type> u_space(L_space * u_raw_space * sigma_space);

  val += GMRF(Q)(u_raw_space);

  vector<Type> mu(beta0 +
                  Z_space * u_space +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

Q <- diag(rowSums(adj.mat)) - adj.mat

Aconstr <- matrix(1, ncol = ncol(Q))
qrc <- qr(t(Aconstr))
L_space <- qr.Q(qrc,complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                L_space = L_space,
                Z_space = Matrix::sparse.model.matrix(~0 + IDf, data),
                Q = as(t(L_space) %*% Q %*% L_space, "dgCMatrix"))

tmbpar <- list(beta0 = 0,
               log_prec_space = 0,
               u_raw_space = numeric(ncol(tmbdata$L_space)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space"),
                      DLL = dll)


tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

tmbfit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

sdr2d <- TMB::sdreport(obj)
summary(sdr2d, "all")


#' Hyper parameters

c(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr2, "fixed")
summary(sdr2d, "fixed")


#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr2, "random")[1, ]
summary(sdr2d, "random")[1, ]


#' Random effects

inlafit$summary.random[[1]]

plot(summary(sdr2, "random")[-1, 1], summary(sdr2d, "report")[ , 1])
abline(0, 1, col = "red")

plot(summary(sdr2, "random")[-1 , 2], summary(sdr2d, "report")[ , 2])
abline(0, 1, col = "red")



#' ## Model 3: RW1 time trend

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_time);
  DATA_SPARSE_MATRIX(R_time); // Structure matrix for RW1

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_time);
  val -= dlgamma(log_prec_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_time(exp(-0.5 * log_prec_time));

  PARAMETER_VECTOR(u_raw_time);
  vector<Type> u_time(u_raw_time * sigma_time);

  val += GMRF(R_time)(u_raw_time);
  val -= dnorm(sum(u_raw_time), Type(0.0), Type(0.001) * u_raw_time.size(), true); // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_time * u_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

Q <- diag(rowSums(adj.mat)) - adj.mat

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), 1e-6)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_time = Matrix::sparse.model.matrix(~0 + Yearf, data),
                R_time = R_time_adj)

tmbpar <- list(beta0 = 0,
               log_prec_time = 0,
               u_raw_time = numeric(ncol(tmbdata$Z_time)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr3 <- TMB::sdreport(obj)
summary(sdr3, "all")

inlafit <- inla(Observed ~ f(Year, model = "rw1", hyper = prec.prior),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))

summary(inlafit)

#' Hyper parameters

c(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr3, "fixed")

#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr3, "random")[1, ]

#' Random effects

inlafit$summary.random[[1]]

plot(inlafit$summary.random[[1]][ , 2], summary(sdr3, "report")[ , 1])
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][ , 3], summary(sdr3, "report")[ , 2])
abline(0, 1, col = "red")


#' #### INLA with custom Cmatrix

inlafitC <- inla(Observed ~ f(ID.Year, model = "generic0", Cmatrix = R_time,
                              hyper = prec.prior, diagonal = diagval, constr = TRUE),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

#' These match:

summary(inlafit)
summary(inlafitC)


#' ### Model 3b: Explicit sum-to-zero constraint

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(L_time); 
  DATA_SPARSE_MATRIX(Z_time);
  DATA_SPARSE_MATRIX(LRL_time); // Structure matrix for RW1

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_time);
  val -= dlgamma(log_prec_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_time(exp(-0.5 * log_prec_time));

  PARAMETER_VECTOR(u_raw_time);
  vector<Type> u_time(L_time * u_raw_time * sigma_time);

  val += GMRF(LRL_time)(u_raw_time);

  vector<Type> mu(beta0 +
                  Z_time * u_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_time);

  return val;
}
'


dll <- tmb_compile_and_load(mod)

Aconstr <- matrix(1, ncol = ncol(R_time))
qrc <- qr(t(Aconstr))
L_time <- qr.Q(qrc, complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                L_time = L_time,
                Z_time = Matrix::sparse.model.matrix(~0 + Yearf, data),
                LRL_time = as(t(L_time) %*% R_time %*% L_time, "dgCMatrix"))

tmbpar <- list(beta0 = 0,
               log_prec_time = 0,
               u_raw_time = numeric(ncol(tmbdata$L_time)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr3b <- TMB::sdreport(obj)
summary(sdr3b, "all")


#' Hyper parameters

c(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr3, "fixed")
summary(sdr3b, "fixed")


#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr3, "random")[1, ]
summary(sdr3b, "random")[1, ]


#' Random effects

inlafit$summary.random[[1]]

plot(summary(sdr3, "report")[ , 1], summary(sdr3b, "report")[ , 1])
abline(0, 1, col = "red")

plot(summary(sdr3, "report")[ , 2], summary(sdr3b, "report")[ , 2])
abline(0, 1, col = "red")



#' ## Model 4: RW1 time + ICAR space

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_time);
  DATA_SPARSE_MATRIX(R_time);

  DATA_SPARSE_MATRIX(Z_space);
  DATA_SPARSE_MATRIX(R_space);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_time);
  val -= dlgamma(log_prec_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_time(exp(-0.5 * log_prec_time));

  PARAMETER_VECTOR(u_raw_time);
  vector<Type> u_time(u_raw_time * sigma_time);

  val += GMRF(R_time)(u_raw_time);
  val -= dnorm(sum(u_raw_time), Type(0.0), Type(0.001) * u_raw_time.size(), true); // soft sum-to-zero constraint

  PARAMETER(log_prec_space);
  val -= dlgamma(log_prec_space, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space(exp(-0.5 * log_prec_space));

  PARAMETER_VECTOR(u_raw_space);
  vector<Type> u_space(u_raw_space * sigma_space);

  val += GMRF(R_space)(u_raw_space);
  val -= dnorm(sum(u_raw_space), Type(0.0), Type(0.001) * u_raw_space.size(), true); // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_time * u_time +
                  Z_space * u_space +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_time);
  ADREPORT(u_space);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

R_space <- diag(rowSums(adj.mat)) - adj.mat
R_space_adj <- R_space + Matrix::Diagonal(ncol(R_space), diagval)

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), diagval)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_time = Matrix::sparse.model.matrix(~0 + Yearf, data),
                R_time = R_time_adj,
                Z_space = Matrix::sparse.model.matrix(~0 + IDf, data),
                R_space = R_space_adj)

tmbpar <- list(beta0 = 0,
               log_prec_time = 0,
               u_raw_time = numeric(ncol(tmbdata$Z_time)),
               log_prec_space = 0,
               u_raw_space = numeric(ncol(tmbdata$Z_space)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_time", "u_raw_space"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr4 <- TMB::sdreport(obj)
summary(sdr4, "all")

inlafit <- inla(Observed ~
                  f(Year, model = "rw1", hyper = prec.prior) +
                  f(as.integer(IDf), model = "besag", hyper = prec.prior, graph = adj.mat),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))


summary(inlafit)

#' Hyper parameters

cbind("mean" = inlafit$misc$theta.mode,
      "se" = sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr4, "fixed")

#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr4, "random")[1, ]

#' Random effects

sdr4sum <- summary(sdr4, "all")

plot(inlafit$summary.random[[1]][ , 2], sdr4sum[rownames(sdr4sum) == "u_time", 1],
     main = "f(Year): mean")
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][ , 3], sdr4sum[rownames(sdr4sum) == "u_time", 2],
     main = "f(Year): sd")
abline(0, 1, col = "red")


plot(inlafit$summary.random[[2]][ , 2], sdr4sum[rownames(sdr4sum) == "u_space", 1],
     main = "f(area): mean")
abline(0, 1, col = "red")

plot(inlafit$summary.random[[2]][ , 3], sdr4sum[rownames(sdr4sum) == "u_space", 2],
     main = "f(area): sd")
abline(0, 1, col = "red")

#' #### INLA with custom Cmatrix

inlafitC <- inla(Observed ~
                   f(ID.Year, model = "generic0", Cmatrix = R_time,
                     hyper = prec.prior, diagonal = diagval, constr = TRUE) +
                   f(as.integer(IDf), model = "generic0", Cmatrix = R_space,
                     hyper = prec.prior, diagonal = diagval, constr = TRUE),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

#' These match:

summary(inlafit)
summary(inlafitC)


#' ### Model 4b: Explicity sum-to-zero constraint

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(L_time);
  DATA_SPARSE_MATRIX(Z_time);
  DATA_SPARSE_MATRIX(LRL_time);

  DATA_MATRIX(L_space);
  DATA_SPARSE_MATRIX(Z_space);
  DATA_SPARSE_MATRIX(LRL_space);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_time);
  val -= dlgamma(log_prec_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_time(exp(-0.5 * log_prec_time));

  PARAMETER_VECTOR(u_raw_time);
  vector<Type> u_time(L_time * u_raw_time * sigma_time);

  val += GMRF(LRL_time)(u_raw_time);

  PARAMETER(log_prec_space);
  val -= dlgamma(log_prec_space, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space(exp(-0.5 * log_prec_space));

  PARAMETER_VECTOR(u_raw_space);
  vector<Type> u_space(L_space * u_raw_space * sigma_space);

  val += GMRF(LRL_space)(u_raw_space);

  vector<Type> mu(beta0 +
                  Z_time * u_time +
                  Z_space * u_space +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_time);
  ADREPORT(u_space);

  return val;
  using namespace density;
}
'

dll <- tmb_compile_and_load(mod)

Aconstr <- matrix(1, ncol = ncol(R_time))
qrc <- qr(t(Aconstr))
L_time <- qr.Q(qrc, complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]

Aconstr <- matrix(1, ncol = ncol(R_space))
qrc <- qr(t(Aconstr))
L_space <- qr.Q(qrc, complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                L_time = L_time,
                Z_time = Matrix::sparse.model.matrix(~0 + Yearf, data),
                LRL_time = as(t(L_time) %*% R_time %*% L_time, "dgCMatrix"),
                L_space = L_space,
                Z_space = Matrix::sparse.model.matrix(~0 + IDf, data),
                LRL_space = as(t(L_space) %*% R_space %*% L_space, "dgCMatrix"))

tmbpar <- list(beta0 = 0,
               log_prec_time = 0,
               u_raw_time = numeric(ncol(tmbdata$L_time)),
               log_prec_space = 0,
               u_raw_space = numeric(ncol(tmbdata$L_space)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_time", "u_raw_space"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr4b <- TMB::sdreport(obj)
summary(sdr4b, "all")


#' Hyper parameters

cbind(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr4, "fixed")
summary(sdr4b, "fixed")


#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr4, "random")[1, ]
summary(sdr4b, "random")[1, ]



#' Random effects

sdr4bsum <- summary(sdr4, "all")

plot(sdr4sum[rownames(sdr4sum) == "u_time", 1],
     sdr4bsum[rownames(sdr4bsum) == "u_time", 1],
     main = "f(Year): mean")
abline(0, 1, col = "red")

plot(sdr4sum[rownames(sdr4sum) == "u_time", 2],
     sdr4bsum[rownames(sdr4bsum) == "u_time", 2],
     main = "f(Year): sd")
abline(0, 1, col = "red")

plot(sdr4sum[rownames(sdr4sum) == "u_space", 1],
     sdr4bsum[rownames(sdr4bsum) == "u_space", 1],
     main = "f(area): mean")
abline(0, 1, col = "red")

plot(sdr4sum[rownames(sdr4sum) == "u_space", 2],
     sdr4bsum[rownames(sdr4bsum) == "u_space", 2],
     main = "f(area): sd")
abline(0, 1, col = "red")



#' ## Model 5: ICAR space x IID time
#'
#' INLA automatically applies a sum-to-zero constraint for the spatial field at each time.
#' This is too many constraints for an intercept plus interaction only model because it
#' imposes that the average is the intercept at each time, and hence there is no time trend
#' in the model.
#'
#' The correct specification for this model should probably be a single sum-to-zero constraint
#' to account for the lost degree of freedom from the intercept term
#' 

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(R_time);
  DATA_SPARSE_MATRIX(R_space);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time * sigma_space_time);

  val += SEPARABLE(GMRF(R_time), GMRF(R_space))(u_raw_space_time);
  for (int i = 0; i < u_raw_space_time.cols(); i++) {
     val -= dnorm(u_raw_space_time.col(i).sum(), Type(0), Type(0.001) * u_raw_space_time.rows(), true);
  }

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'


dll <- tmb_compile_and_load(mod)


R_space <- diag(rowSums(adj.mat)) - adj.mat
R_space_scaled <- inla.scale.model(R_space, constr = list(A = matrix(1, ncol = ncol(R_space)), e = 0))
R_space_adj <- R_space + Matrix::Diagonal(ncol(R_space), 1e-6)
R_space_scaled_adj <- R_space_scaled + Matrix::Diagonal(ncol(R_space_scaled), 1e-6)

R_time <- Matrix::sparseMatrix(1:19, 1:19, x = rep(1L, 19))

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_time = R_time,
                R_space = R_space_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr5 <- TMB::sdreport(obj)
summary(sdr5, "all")

inlafit <- inla(Observed ~
                  f(as.integer(IDf), model = "besag", hyper = prec.prior, scale.model = FALSE,
                    graph = adj.mat, diagonal = 1e-6,
                    group = ID.Year, control.group = list(model = "iid")),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))

summary(inlafit)


#' Hyper parameter comparison
#'

cbind("mean" = inlafit$misc$theta.mode,
      "se" = sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr5, "fixed")


#' Fixed effects (Intercept)

inlafit$summary.fixed[ , 1:2]
summary(sdr5, "random")[1, , drop = FALSE]


#' Random effects mean and standard deviation

plot(inlafit$summary.random[[1]][,2], summary(sdr5, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][,3], summary(sdr5, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1, col = "red")


#' #### INLA with custom Cmatrix

R_space_time <- kronecker(R_time, R_space)

Aconstr <- t(model.matrix(~0+factor(ID.Year), data[order(data$id.area.year), ]))

inlafitC <- inla(Observed ~
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time,
                     hyper = prec.prior, diagonal = 1e-6,
                     extraconstr = list(A = Aconstr, e = numeric(nrow(Aconstr)))),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))


#' These match:

summary(inlafit)
summary(inlafitC)

#' ### Model 5b: Express GMRF as kronecker product

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(Aconstr)
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(R_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_raw_space_time_v(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time_v * sigma_space_time);

  val += GMRF(R_space_time)(u_raw_space_time);
  val -= dnorm(Aconstr * u_raw_space_time_v, Type(0), Type(0.001) * u_raw_space_time.rows(), true).sum();  // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

R_space_time_adj <- R_space_time + Matrix::Diagonal(ncol(R_space_time), diagval)


tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Aconstr = Aconstr,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_space_time = R_space_time_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr5b <- TMB::sdreport(obj)
summary(sdr5b, "all")

#' Hyper parmaters 

summary(sdr5, "fixed")
summary(sdr5b, "fixed")

#' Intercept

summary(sdr5, "random")[1, ]
summary(sdr5b, "random")[1, ]


#' Random effects

plot(summary(sdr5, "report")[ , 1],
     summary(sdr5b, "report")[ , 1])
abline(0, 1, col = "red")

plot(summary(sdr5, "report")[ , 2],
     summary(sdr5b, "report")[ , 2])
abline(0, 1, col = "red")

#' ### Model 5c: explicit sum-to-zero constraint

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(L_space_time)
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(LRL_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_VECTOR(u_raw_space_time);
  vector<Type> u_space_time(L_space_time * u_raw_space_time * sigma_space_time);

  val += GMRF(LRL_space_time)(u_raw_space_time);

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

qrc <- qr(t(Aconstr))
L_space_time <- qr.Q(qrc, complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]


R_space_time_adj <- R_space_time + Matrix::Diagonal(ncol(R_space_time), diagval)


tmbdata <- list(y = data$Observed,
                E = data$Expected,
                L_space_time = L_space_time,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                LRL_space_time = as(t(L_space_time) %*% R_space_time %*% L_space_time, "dgCMatrix"))

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = numeric(ncol(tmbdata$L_space_time)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr5c <- TMB::sdreport(obj)
summary(sdr5c, "all")



#' Hyper parameters

cbind(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr5, "fixed")
summary(sdr5c, "fixed")


#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr5, "random")[1, ]
summary(sdr5c, "random")[1, ]



#' Random effects

sdr5csum <- summary(sdr5, "all")

plot(summary(sdr5, "report")[ , 1],
     summary(sdr5c, "report")[ , 1])
abline(0, 1, col = "red")

plot(summary(sdr5, "report")[ , 2],
     summary(sdr5c, "report")[ , 2])
abline(0, 1, col = "red")



#' ## Model 6: IID x RW1
#'
#' Use and IID spatial effect, but a RW1 time effect. In this case, INLA by default does not
#' add any constraint because the main effect is proper.
#'
#' To match the `Cmatrix` implementation with the IID x RW1 implementation with the
#' `group` option, the argument `rankdef` needs to be specified as the rank deficiency
#' of the Kronecker product precision matrix. This is the number of areas.
#' 

inlafit6 <- inla(Observed ~ 
                   f(as.integer(IDf), model = "iid", hyper = prec.prior,
                     graph = adj.mat, diagonal = diagval,  
                     group = ID.Year, control.group = list(model = "rw1", scale.model = TRUE)),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

#' Note that rank deficiency line in the log file pertains to the rank deficiency of the
#' main IID effect, which has zero rank deficiency:
#'

grep(".*rank.*", inlafit6$logfile, value = TRUE)

#' For the `group` model, the INLA log file does not make any comment about the rank
#' deficiency. However, we will see below that internally that the rank deficiency
#' which is imposed RW1 group model is accounted for.

#' #### INLA with custom Cmatrix


R_space <- diag(ncol(adj.mat))

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), 1e-6)

R_time_scaled <- inla.scale.model(R_time, constr = list(A = matrix(1, ncol = ncol(R_time)), e = 0))

R_space_time <- kronecker(R_time_scaled, R_space)

inlafit6C <- inla(Observed ~
                    f(id.area.year, model = "generic0", Cmatrix = R_space_time, hyper = prec.prior,
                      diagonal = diagval, rankdef = 32),
                  data = data, E = Expected, family = "poisson",
                  control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                  control.compute = list(config = TRUE))

#' Here's the `rankdef` comment in the logfile:

grep(".*rank.*", inlafit6C$logfile, value = TRUE)

#' These results match:

summary(inlafit6)
summary(inlafit6C)


#' Marginal likelihood

inlafit6$mlik
inlafit6C$mlik


#' Same constraints (none)

all.equal(inlafit6$misc$configs$constr,
          inlafit6C$misc$configs$constr)


#' Hyper parameters

inlafit6$internal.summary.hyperpar[, 1:2]
inlafit6C$internal.summary.hyperpar[, 1:2]

#' Fixed effects

inlafit6$summary.fixed
inlafit6C$summary.fixed


#' Random effects

plot(inlafit6$summary.random[[1]][,2],
     inlafit6C$summary.random[[1]][,2],
     main = "Random effect mean")
abline(a = 0, b = 1, col = "red")

plot(inlafit6$summary.random[[1]][,3],
     inlafit6C$summary.random[[1]][,3],
     main = "Random effect standard deviation")
abline(a = 0, b = 1, col = "red")



#' ### TMB 

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(R_time);
  DATA_SPARSE_MATRIX(R_space);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time * sigma_space_time);

  val += SEPARABLE(GMRF(R_time), GMRF(R_space))(u_raw_space_time);

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_space = as(R_space, "dgTMatrix"),
                R_time = R_time_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr6 <- TMB::sdreport(obj)
summary(sdr6, "all")


#' Hyper parameter comparison
#'
#' The standard devation from TMB is somewhat larger than INLA.

summary(sdr6, "fixed")

cbind("mean" = inlafit6$misc$theta.mode,
      "se" = sqrt(diag(inlafit6$misc$cov.intern)))


#' Fixed effects (Intercept)

summary(sdr6, "random")[1, , drop = FALSE]
inlafit6$summary.fixed[ , 1:2]
inlafit6C$summary.fixed[ , 1:2]


#' Random effects mean and standard deviation

plot(inlafit6$summary.random[[1]][,2], summary(sdr6, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1)

plot(inlafit6$summary.random[[1]][,3], summary(sdr6, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1)



#' ### Model 6b: Express GMRF as a Kronecker product

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(R_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time * sigma_space_time);

  val += GMRF(R_space_time)(u_raw_space_time);

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

## R_space_time_adj <- R_space_time + Matrix::Diagonal(ncol(R_space_time), diagval)
R_space_time_adj <- kronecker(R_time_adj, R_space)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_space_time = R_space_time_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr6b <- TMB::sdreport(obj)
summary(sdr6b, "all")


#' Hyper parameter comparison
#'
#' The standard devation from TMB is somewhat larger than INLA.

summary(sdr6, "fixed")
summary(sdr6b, "fixed")

cbind("mean" = inlafit6$misc$theta.mode,
      "se" = sqrt(diag(inlafit6$misc$cov.intern)))


#' Fixed effects (Intercept)

summary(sdr6, "random")[1, , drop = FALSE]
summary(sdr6b, "random")[1, , drop = FALSE]


#' Random effects mean and standard deviation

plot(summary(sdr6, "report")[,1], summary(sdr6b, "report")[,1],
     xlab = "TMB SEPARABLE()", ylab = "TMB Kronecker()", main = "Random effect point estimates")
abline(0, 1)

plot(summary(sdr6, "report")[,2], summary(sdr6b, "report")[,2],
     xlab = "TMB SEPARABLE()", ylab = "TMB Kronecker()", main = "Random effect standard devation")
abline(0, 1)


#' ## Model 7: IID x RW1 with constraints
#'
#' By changing the order of the main and group terms in the INLA smooth formula,
#' we can implement the same model with constratints on the RW1 terms
#'
#' First fit the unconstrained model and confirm it matches the previous version.

inlafit7_unconstr <- inla(Observed ~ 
                            f(ID.Year, model = "rw1", hyper = prec.prior,
                              scale.model = TRUE, diagonal = diagval, constr = FALSE,
                              group = ID, control.group = list(model = "iid")),
                          data = data, E = Expected, family = "poisson",
                          control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                          control.compute = list(config = TRUE))

summary(inlafit6)
summary(inlafit7_unconstr)

#' Fit the model with default sum-to-zero constraints on each group.

inlafit7 <- inla(Observed ~ 
                   f(ID.Year, model = "rw1", hyper = prec.prior,
                     scale.model = TRUE, diagonal = diagval, constr = TRUE,
                     group = ID, control.group = list(model = "iid")),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

grep(".*rank.*", inlafit7_unconstr$logfile, value = TRUE)


#' #### INLA with custom Cmatrix


R_space <- diag(ncol(adj.mat))

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), 1e-6)

R_time_scaled <- inla.scale.model(R_time, constr = list(A = matrix(1, ncol = ncol(R_time)), e = 0))

R_space_time <- kronecker(R_time_scaled, R_space)

Aconstr <- t(model.matrix(~0+factor(IDf), data[order(data$id.area.year), ]))

inlafit7C <- inla(Observed ~
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time, hyper = prec.prior,
                     extraconstr = list(A = Aconstr, e = numeric(nrow(Aconstr))),
                     diagonal = diagval),
                  data = data, E = Expected, family = "poisson",
                  control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                  control.compute = list(config = TRUE))


grep(".*rank.*", inlafit7C$logfile, value = TRUE)

#' These results match:

summary(inlafit7)
summary(inlafit7C)


#' Marginal likelihood

inlafit7$mlik
inlafit7C$mlik


#' Same constraints (none)

dim(inlafit7$misc$configs$constr$A)
dim(inlafit7C$misc$configs$constr$A)


#' Hyper parameters

inlafit7$internal.summary.hyperpar[, 1:2]
inlafit7C$internal.summary.hyperpar[, 1:2]

#' Fixed effects

inlafit7$summary.fixed
inlafit7C$summary.fixed


#' Random effects

plot(inlafit7$summary.random[[1]][order(inlafit7$summary.random[[1]]$ID), 2],
     inlafit7C$summary.random[[1]][,2],
     main = "Random effect mean")
abline(a = 0, b = 1, col = "red")

plot(inlafit7$summary.random[[1]][order(inlafit7$summary.random[[1]]$ID), 3],
     inlafit7C$summary.random[[1]][,3],
     main = "Random effect standard deviation")
abline(a = 0, b = 1, col = "red")



#' ### TMB 

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(Aconstr);
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(R_time);
  DATA_SPARSE_MATRIX(R_space);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_raw_space_time_v(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time_v * sigma_space_time);

  val += SEPARABLE(GMRF(R_time), GMRF(R_space))(u_raw_space_time);
  val -= dnorm(Aconstr * u_raw_space_time_v, Type(0), Type(0.001) * u_raw_space_time.cols(), true).sum();  // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Aconstr = Aconstr,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_space = as(R_space, "dgTMatrix"),
                R_time = R_time_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr7 <- TMB::sdreport(obj)
summary(sdr7, "all")


#' Hyper parameter comparison
#'
#' The standard devation from TMB is somewhat larger than INLA.

summary(sdr7, "fixed")

cbind("mean" = inlafit7C$misc$theta.mode,
      "se" = sqrt(diag(inlafit7C$misc$cov.intern)))


#' Fixed effects (Intercept)

summary(sdr7, "random")[1, , drop = FALSE]
inlafit7$summary.fixed[ , 1:2]
inlafit7C$summary.fixed[ , 1:2]


#' Random effects mean and standard deviation

plot(inlafit7C$summary.random[[1]][,2], summary(sdr7, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(inlafit7C$summary.random[[1]][,3], summary(sdr7, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1, col = "red")



#' ### Model 7b: Explicit sum-to-zero constraint

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(L_space_time)
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(LRL_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_VECTOR(u_raw_space_time);
  vector<Type> u_space_time(L_space_time * u_raw_space_time * sigma_space_time);

  val += GMRF(LRL_space_time)(u_raw_space_time);

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

qrc <- qr(t(Aconstr))
L_space_time <- qr.Q(qrc, complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                L_space_time = L_space_time,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                LRL_space_time = as(t(L_space_time) %*% R_space_time %*% L_space_time, "dgCMatrix"))

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = numeric(ncol(tmbdata$L_space_time)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)


sdr7b <- TMB::sdreport(obj)
summary(sdr7b, "all")


#' Hyper parameter comparison


summary(sdr7, "fixed")
summary(sdr7b, "fixed")

cbind("mean" = inlafit7C$misc$theta.mode,
      "se" = sqrt(diag(inlafit7C$misc$cov.intern)))


#' Fixed effects (Intercept)

summary(sdr7, "random")[1, , drop = FALSE]
summary(sdr7b, "random")[1, , drop = FALSE]

inlafit7C$summary.fixed[ , 1:2]

#' Random effects mean and standard deviation

plot(summary(sdr7, "report")[,1], summary(sdr7b, "report")[,1],
     xlab = "TMB soft constraint)", ylab = "TMB hard constraint",
     main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(summary(sdr7, "report")[,2], summary(sdr7b, "report")[,2],
     xlab = "TMB soft constraint", ylab = "TMB hard constraint",
     main = "Random effect standard devation")
abline(0, 1, col = "red")

plot(inlafit7C$summary.random[[1]][ , 2], summary(sdr7b, "report")[,1],
     xlab = "TMB soft constraint)", ylab = "TMB hard constraint",
     main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(inlafit7C$summary.random[[1]][ , 3], summary(sdr7b, "report")[,2],
     xlab = "TMB soft constraint", ylab = "TMB hard constraint",
     main = "Random effect standard devation")
abline(0, 1, col = "red")



#' ## Model 8: ICAR space x RW1
#'
#' Change the time trend to a RW1. The TMB code remains the same. The only change is the
#' structure matrix for the time trend. If we were using a scaled parameterisation, we
#' would need to adjust the penalty for the additional rank deficiency of the structure
#' matrix.
#' 
#' For the `Cmatrix` version of the of the INLA model, the rank deficiency is $rows +
#' columns - 1$ since both are rank deficient one.
#'
#' INLA automatically applies a sum-to-zero constraint for the spatial field at each time.
#' This is too many constraints for an intercept plus interaction only model because it
#' imposes that the average is the intercept at each time, and hence there is no time trend
#' in the model.
#'
#' The correct specification for this model should probably be a single sum-to-zero constraint
#' to account for the lost degree of freedom from the intercept term
#' 

inlafit8 <- inla(Observed ~ 
                   f(Year, model = "rw1", scale.model = TRUE,
                     hyper = prec.prior, diagonal = diagval, constr = TRUE,
                     group = ID,
                     control.group = list(model = "besag", graph = adj.mat, scale.model = TRUE)),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

#' Only the main effect (ICAR) rank deficiency is reflected in the log file.

grep(".*rank.*", inlafit8$logfile, value = TRUE)

#' #### INLA with custom Cmatrix


R_space <- diag(rowSums(adj.mat)) - adj.mat
R_space_scaled <- inla.scale.model(R_space, constr = list(A = matrix(1, ncol = ncol(R_space)), e = 0))
R_space_scaled_adj <- R_space_scaled + Matrix::Diagonal(ncol(R_space_scaled), diagval)

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_scaled <- inla.scale.model(R_time, constr = list(A = matrix(1, ncol = ncol(R_time)), e = 0))
R_time_scaled_adj <- R_time_scaled + Matrix::Diagonal(ncol(R_time_scaled), diagval)

R_space_time <- kronecker(R_time_scaled, R_space_scaled)

rankdef <- nrow(R_space) + ncol(R_time) - 1

Aconstr <- t(model.matrix(~0+IDf, data[order(data$id.area.year), ]))

inlafit8C <- inla(Observed ~
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time,
                     hyper = prec.prior, diagonal = diagval, rankdef = rankdef, 
                     extraconstr = list(A = Aconstr, e = numeric(nrow(Aconstr)))),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

grep(".*rank.*", inlafit8C$logfile, value = TRUE)


#' The precision parameter estimates match, but the marginal likelihood does not match.

summary(inlafit8)
summary(inlafit8C)

#' Marginal likelihood

inlafit8$mlik
inlafit8C$mlik


#' Same constraints (none)

all.equal(inlafit8$misc$configs$constr,
          inlafit8C$misc$configs$constr)


#' Hyper parameters

inlafit8$internal.summary.hyperpar[, 1:2]
inlafit8C$internal.summary.hyperpar[, 1:2]


#' Fixed effects

inlafit8$summary.fixed[ , 1:2]
inlafit8C$summary.fixed[ , 1:2]

#' Random effects

plot(inlafit8$summary.random[[1]][order(inlafit8$summary.random[[1]]$ID), 2],
     inlafit8C$summary.random[[1]][,2],
     main = "Random effect mean")
abline(a = 0, b = 1, col = "red")

plot(inlafit8$summary.random[[1]][order(inlafit8$summary.random[[1]]$ID), 3],
     inlafit8C$summary.random[[1]][,3],
     main = "Random effect standard deviation")
abline(a = 0, b = 1, col = "red")

#' ### TMB 

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(Aconstr);
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(R_time);
  DATA_SPARSE_MATRIX(R_space);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_raw_space_time_v(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time * sigma_space_time);

  val += SEPARABLE(GMRF(R_time), GMRF(R_space))(u_raw_space_time);
  val -= dnorm(Aconstr * u_raw_space_time_v, Type(0), Type(0.001) * u_raw_space_time.cols(), true).sum();  // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

R_space_time_adj <- R_space_time + Matrix::Diagonal(ncol(R_space_time), 1e-6)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Aconstr = Aconstr,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_space = R_space_scaled_adj,
                R_time = R_time_scaled_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr8 <- TMB::sdreport(obj)
summary(sdr8, "all")


#' Hyper parameter comparison


summary(sdr8, "fixed")

cbind("mean" = inlafit8$misc$theta.mode,
      "se" = sqrt(diag(inlafit8$misc$cov.intern)))

cbind("mean" = inlafit8C$misc$theta.mode,
      "se" = sqrt(diag(inlafit8C$misc$cov.intern)))



#' Fixed effects (Intercept)

summary(sdr8, "random")[1, , drop = FALSE]
inlafit8$summary.fixed[ , 1:2]
inlafit8C$summary.fixed[ , 1:2]



#' Random effects mean and standard deviation

plot(inlafit8C$summary.random[[1]][,2], summary(sdr8, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(inlafit8C$summary.random[[1]][,3], summary(sdr8, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1, col = "red")


#' ### Model 8b: Explicit sum-to-zero constraint

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(L_space_time)
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(LRL_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_VECTOR(u_raw_space_time);
  vector<Type> u_space_time(L_space_time * u_raw_space_time * sigma_space_time);

  val += GMRF(LRL_space_time)(u_raw_space_time);

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

qrc <- qr(t(Aconstr))
L_space_time <- qr.Q(qrc, complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]

LRL <- as(t(L_space_time) %*% R_space_time %*% L_space_time, "dgCMatrix")
LRL_adj <- LRL + Matrix::Diagonal(ncol(LRL), diagval)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                L_space_time = L_space_time,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                LRL_space_time = LRL_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = numeric(ncol(tmbdata$L_space_time)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)


sdr8b <- TMB::sdreport(obj)
summary(sdr8b, "all")


#' Hyper parameter comparison


summary(sdr8, "fixed")
summary(sdr8b, "fixed")

cbind("mean" = inlafit8C$misc$theta.mode,
      "se" = sqrt(diag(inlafit8C$misc$cov.intern)))


#' Fixed effects (Intercept)

summary(sdr8, "random")[1, , drop = FALSE]
summary(sdr8b, "random")[1, , drop = FALSE]

inlafit8C$summary.fixed[ , 1:2]

#' Random effects mean and standard deviation

plot(summary(sdr8, "report")[,1], summary(sdr8b, "report")[,1],
     xlab = "TMB soft constraint)", ylab = "TMB hard constraint",
     main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(summary(sdr8, "report")[,2], summary(sdr8b, "report")[,2],
     xlab = "TMB soft constraint", ylab = "TMB hard constraint",
     main = "Random effect standard devation")
abline(0, 1, col = "red")

plot(inlafit8C$summary.random[[1]][ , 2], summary(sdr8b, "report")[,1],
     xlab = "TMB soft constraint)", ylab = "TMB hard constraint",
     main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(inlafit8C$summary.random[[1]][ , 3], summary(sdr8b, "report")[,2],
     xlab = "TMB soft constraint", ylab = "TMB hard constraint",
     main = "Random effect standard devation")
abline(0, 1, col = "red")


#' ## Model 9: ICAR space x RW1 time with full constraints
#'


R_space <- diag(rowSums(adj.mat)) - adj.mat
R_space_scaled <- inla.scale.model(R_space, constr = list(A = matrix(1, ncol = ncol(R_space)), e = 0))
R_space_scaled_adj <- R_space_scaled + Matrix::Diagonal(ncol(R_space_scaled), diagval)

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_scaled <- inla.scale.model(R_time, constr = list(A = matrix(1, ncol = ncol(R_time)), e = 0))
R_time_scaled_adj <- R_time_scaled + Matrix::Diagonal(ncol(R_time_scaled), diagval)

R_space_time <- kronecker(R_time_scaled, R_space_scaled)

eig <- eigen(R_space_time)
eigval <- zapsmall(eig$values)
eigvec <- eig$vectors

rankdef <- nrow(R_space) + ncol(R_time) - 1

sum(eigval == 0)
rankdef

Aconstr <- t(eigvec[ , eigval == 0])

inlafit9C <- inla(Observed ~
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time,
                     hyper = prec.prior, diagonal = diagval, rankdef = rankdef, 
                     extraconstr = list(A = Aconstr, e = numeric(nrow(Aconstr)))),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

grep(".*rank.*", inlafit9C$logfile, value = TRUE)


#' The precision parameter estimates match, but the marginal likelihood does not match.

summary(inlafit9C)

#' Marginal likelihood

inlafit9$mlik
inlafit9C$mlik


#' Same constraints (none)

all.equal(inlafit9$misc$configs$constr,
          inlafit9C$misc$configs$constr)


#' Hyper parameters

inlafit9$internal.summary.hyperpar[, 1:2]
inlafit9C$internal.summary.hyperpar[, 1:2]


#' Fixed effects

inlafit9$summary.fixed[ , 1:2]
inlafit9C$summary.fixed[ , 1:2]

#' Random effects

plot(inlafit9$summary.random[[1]][order(inlafit9$summary.random[[1]]$ID), 2],
     inlafit9C$summary.random[[1]][,2],
     main = "Random effect mean")
abline(a = 0, b = 1, col = "red")

plot(inlafit9$summary.random[[1]][order(inlafit9$summary.random[[1]]$ID), 3],
     inlafit9C$summary.random[[1]][,3],
     main = "Random effect standard deviation")
abline(a = 0, b = 1, col = "red")

#' ### TMB 

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(Aconstr);
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(R_time);
  DATA_SPARSE_MATRIX(R_space);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_raw_space_time_v(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time * sigma_space_time);

  val += SEPARABLE(GMRF(R_time), GMRF(R_space))(u_raw_space_time);
  val -= dnorm(Aconstr * u_raw_space_time_v, Type(0), Type(0.001) * u_raw_space_time.cols(), true).sum();  // soft sum-to-zero constraint

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

R_space_time_adj <- R_space_time + Matrix::Diagonal(ncol(R_space_time), 1e-6)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Aconstr = Aconstr,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_space = R_space_scaled_adj,
                R_time = R_time_scaled_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)

sdr9 <- TMB::sdreport(obj)
summary(sdr9, "all")


#' Hyper parameter comparison


summary(sdr9, "fixed")

cbind("mean" = inlafit9C$misc$theta.mode,
      "se" = sqrt(diag(inlafit9C$misc$cov.intern)))



#' Fixed effects (Intercept)

summary(sdr9, "random")[1, , drop = FALSE]
inlafit9C$summary.fixed[ , 1:2]


#' Random effects mean and standard deviation

plot(inlafit9C$summary.random[[1]][,2], summary(sdr9, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(inlafit9C$summary.random[[1]][,3], summary(sdr9, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1, col = "red")


#' ### Model 9b: Explicit sum-to-zero constraint

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_MATRIX(L_space_time)
  DATA_SPARSE_MATRIX(Z_space_time);
  DATA_SPARSE_MATRIX(LRL_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_VECTOR(u_raw_space_time);
  vector<Type> u_space_time(L_space_time * u_raw_space_time * sigma_space_time);

  val += GMRF(LRL_space_time)(u_raw_space_time);

  vector<Type> mu(beta0 +
                  Z_space_time * u_space_time +
                  log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(u_space_time);

  return val;
}
'

dll <- tmb_compile_and_load(mod)

qrc <- qr(t(Aconstr))
L_space_time <- qr.Q(qrc, complete=TRUE)[ , (nrow(Aconstr)+1):ncol(Aconstr)]

LRL <- as(t(L_space_time) %*% R_space_time %*% L_space_time, "dgCMatrix")
LRL_adj <- LRL + Matrix::Diagonal(ncol(LRL), diagval)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                L_space_time = L_space_time,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                LRL_space_time = LRL_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = numeric(ncol(tmbdata$L_space_time)))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr)


sdr9b <- TMB::sdreport(obj)
summary(sdr9b, "all")


#' Hyper parameter comparison


summary(sdr9, "fixed")
summary(sdr9b, "fixed")

cbind("mean" = inlafit9C$misc$theta.mode,
      "se" = sqrt(diag(inlafit9C$misc$cov.intern)))


#' Fixed effects (Intercept)

summary(sdr9, "random")[1, , drop = FALSE]
summary(sdr9b, "random")[1, , drop = FALSE]

inlafit9C$summary.fixed[ , 1:2]

#' Random effects mean and standard deviation

plot(summary(sdr9, "report")[,1], summary(sdr9b, "report")[,1],
     xlab = "TMB soft constraint)", ylab = "TMB hard constraint",
     main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(summary(sdr9, "report")[,2], summary(sdr9b, "report")[,2],
     xlab = "TMB soft constraint", ylab = "TMB hard constraint",
     main = "Random effect standard devation")
abline(0, 1, col = "red")

