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

data$IDf <- factor(data$ID, 1:max(as.integer(data$ID)))
data$Yearf <- factor(data$Year, unique(data$Year))

data$ID.Year <- data$Year - 1973 + 1
data$ID2 <- data$ID


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


TMB::sdreport(obj)

inlafit <- inla(Observed ~ 1,
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"))

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


#' ### Hyper parameter comparison
#'

cbind("mean" = inlafit$misc$theta.mode,
      "se" = sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr5, "fixed")


#' ### Fixed effects (Intercept)

inlafit$summary.fixed[ , 1:2]
summary(sdr5, "random")[1, , drop = FALSE]


#' ### Random effects mean and standard deviation

plot(inlafit$summary.random[[1]][,2], summary(sdr5, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][,3], summary(sdr5, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1, col = "red")


#' #### INLA with custom Cmatrix

data$area.year <- interaction(data$IDf, data$Year)
data$id.area.year <- as.integer(data$area.year)

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


#' ## Model 6: IID x RW1
#'
#' Use and IID spatial effect, but a RW1 time effect. In this case, INLA by default does not
#' add any constraint because the main effect is proper.
#'
#' To match the `Cmatrix` implementation with the IID x RW1 implementation with the
#' `group` option, the argument `rankdef` needs to be specified as the rank deficiency
#' of the Kronecker product precision matrix. This is the number of areas.
#' 

inlafit <- inla(Observed ~ 
                  f(as.integer(IDf), model = "iid", hyper = prec.prior,
                    graph = adj.mat, diagonal = diagval,  
                    group = ID.Year, control.group = list(model = "rw1", scale.model = TRUE)),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))

#' Note that rank deficiency line in the log file pertains to the rank deficiency of the
#' main IID effect, which has zero rank deficiency:
#'

grep(".*rank.*", inlafit$logfile, value = TRUE)

#' For the `group` model, the INLA log file does not make any comment about the rank
#' deficiency. However, we will see below that internally that the rank deficiency
#' which is imposed RW1 group model is accounted for.

#' #### INLA with custom Cmatrix

data$area.year <- interaction(data$IDf, data$Year)
data$id.area.year <- as.integer(data$area.year)

R_space <- diag(ncol(adj.mat))

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), 1e-6)

R_time_scaled <- inla.scale.model(R_time, constr = list(A = matrix(1, ncol = ncol(R_time)), e = 0))

R_space_time <- kronecker(R_time_scaled, R_space)

inlafitC <- inla(Observed ~
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time, hyper = prec.prior,
                     diagonal = diagval, rankdef = 32),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

#' Here's the `rankdef` comment in the logfile:

grep(".*rank.*", inlafitC$logfile, value = TRUE)

#' These results match:

summary(inlafit)
summary(inlafitC)


#' Marginal likelihood

inlafit$mlik
inlafitC$mlik


#' Same constraints (none)

all.equal(inlafit$misc$configs$constr,
          inlafitC$misc$configs$constr)


#' Hyper parameters

inlafit$internal.summary.hyperpar[, 1:2]
inlafitC$internal.summary.hyperpar[, 1:2]

#' Fixed effects

inlafit$summary.fixed
inlafitC$summary.fixed


#' Random effects

plot(inlafit$summary.random[[1]][,2],
     inlafitC$summary.random[[1]][,2],
     main = "Random effect mean")

plot(inlafit$summary.random[[1]][,3],
     inlafitC$summary.random[[1]][,3],
     main = "Random effect standard deviation")
abline(a = 0, b = 1, col = "red")



#' #### TMB 

mod <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Z_space_time);
  // DATA_SPARSE_MATRIX(R_time);
  // DATA_SPARSE_MATRIX(R_space);
  DATA_SPARSE_MATRIX(R_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time * sigma_space_time);

  // val += SEPARABLE(GMRF(R_time), GMRF(R_space))(u_raw_space_time);
  val += GMRF(R_space_time)(u_raw_space_time);
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

R_space_time_adj <- R_space_time + Matrix::Diagonal(ncol(R_space_time), 1e-6)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                ## R_space = R_space_adj,
                ## R_time = R_time_adj,
                R_space_time = R_space_time_adj)

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


#' ### Hyper parameter comparison
#'
#' The standard devation from TMB is somewhat larger than INLA.

summary(sdr6, "fixed")

cbind("mean" = inlafit$misc$theta.mode,
      "se" = sqrt(diag(inlafit$misc$cov.intern)))


#' ### Fixed effects (Intercept)

summary(sdr6, "random")[1, , drop = FALSE]
inlafit$summary.fixed[ , 1:2]
inlafitC$summary.fixed[ , 1:2]


#' ### Random effects mean and standard deviation

plot(inlafit$summary.random[[1]][,2], summary(sdr6, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1)

plot(inlafit$summary.random[[1]][,3], summary(sdr6, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1)






#' ## Model 7: ICAR space x RW1
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

inlafit <- inla(Observed ~ 
                  f(as.integer(IDf), model = "besag", hyper = prec.prior,
                    graph = adj.mat, diagonal = 1e-6, scale.model = FALSE, constr = TRUE,
                    group = ID.Year, control.group = list(model = "rw1", scale.model = FALSE)),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))

#' Only the main effect (ICAR) rank deficiency is reflected in the log file.

grep(".*rank.*", inlafit$logfile, value = TRUE)

#' #### INLA with custom Cmatrix

data$area.year <- interaction(data$IDf, data$Year)
data$id.area.year <- as.integer(data$area.year)

R_space <- diag(rowSums(adj.mat)) - adj.mat

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), 1e-6)

R_space_time <- kronecker(R_time, R_space)

rankdef <- nrow(R_space) + ncol(R_time) - 1

## Aconstr <- t(model.matrix(~0+factor(ID.Year), data[order(data$id.area.year), ]))
Aconstr <- t(apply(cbind(0:18, 1, 18:0)*32, 1, rep, x = c(0, 1, 0)))

inlafitC <- inla(Observed ~
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time,
                     hyper = prec.prior, diagonal = 1e-6, rankdef = rankdef, 
                     extraconstr = list(A = Aconstr, e = numeric(nrow(Aconstr)))),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))

grep(".*rank.*", inlafitC$logfile, value = TRUE)


#' The precision parameter estimates match, but the marginal likelihood does not match.

summary(inlafit)
summary(inlafitC)

#' Marginal likelihood

inlafit$mlik
inlafitC$mlik


#' Same constraints (none)

all.equal(inlafit$misc$configs$constr,
          inlafitC$misc$configs$constr)


#' Hyper parameters

inlafit$internal.summary.hyperpar[, 1:2]
inlafitC$internal.summary.hyperpar[, 1:2]


#' Fixed effects

inlafit$summary.fixed[ , 1:2]
inlafitC$summary.fixed[ , 1:2]

#' Random effects

plot(inlafit$summary.random[[1]][,2],
     inlafitC$summary.random[[1]][,2],
     main = "Random effect mean")
abline(a = 0, b = 1, col = "red")

plot(inlafit$summary.random[[1]][,3],
     inlafitC$summary.random[[1]][,3],
      xlim = c(2.2, 2.6),
      ylim = c(2.2, 2.6),
     main = "Random effect standard deviation")
abline(a = 0, b = 1, col = "red")

#' #### TMB 

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
  // DATA_SPARSE_MATRIX(R_space_time);

  Type val(0);

  PARAMETER(beta0);
  // beta0 ~ 1

  PARAMETER(log_prec_space_time);
  val -= dlgamma(log_prec_space_time, Type(0.001), Type(1.0 / 0.001), true);
  Type sigma_space_time(exp(-0.5 * log_prec_space_time));

  PARAMETER_ARRAY(u_raw_space_time);
  vector<Type> u_space_time(u_raw_space_time * sigma_space_time);

  val += SEPARABLE(GMRF(R_time), GMRF(R_space))(u_raw_space_time);
  // val += GMRF(R_space_time)(u_raw_space_time);
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

R_space_time_adj <- R_space_time + Matrix::Diagonal(ncol(R_space_time), 1e-6)

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_space = R_space_adj,
                R_time = R_time_adj)

tmbpar <- list(beta0 = 0,
               log_prec_space_time = 0,
               u_raw_space_time = array(0, c(nrow(tmbdata$R_space), nrow(tmbdata$R_time))))

obj <- TMB::MakeADFun(data = tmbdata,
                      parameters = tmbpar,
                      random = c("beta0", "u_raw_space_time"),
                      DLL = dll)

tmbfit <- nlminb(obj$par, obj$fn, obj$gr,
                 control = list(iter.max = 1000,
                                eval.max = 1000))

sdr7 <- TMB::sdreport(obj)
summary(sdr7, "all")


#' ### Hyper parameter comparison
#'
#' The standard devation from TMB is somewhat larger than INLA.

summary(sdr7, "fixed")

cbind("mean" = inlafit$misc$theta.mode,
      "se" = sqrt(diag(inlafit$misc$cov.intern)))

cbind("mean" = inlafitC$misc$theta.mode,
      "se" = sqrt(diag(inlafitC$misc$cov.intern)))



#' ### Fixed effects (Intercept)

summary(sdr7, "random")[1, , drop = FALSE]
inlafit$summary.fixed[ , 1:2]
inlafitC$summary.fixed[ , 1:2]



#' ### Random effects mean and standard deviation

plot(inlafit$summary.random[[1]][,2], summary(sdr7, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1)

plot(inlafit$summary.random[[1]][,3], summary(sdr7, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1)

