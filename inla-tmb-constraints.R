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

data$IDf <- factor(data$ID, 1:max(as.integer(data$ID)))

data2 <- filter(data, ID %in% 1:3, Year %in% 1973:1974) %>%
  droplevels()
adj.mat2 <- Matrix::spMatrix(3, 3, c(1, 2, 2, 3), c(2, 1, 3, 2), rep(1L, 4))
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

sdr <- TMB::sdreport(obj)
summary(sdr, "all")

prec.prior <- list(prec = list(param = c(0.001, 0.001)))

inlafit <- inla(Observed ~ f(as.integer(IDf), model = "besag",
                             hyper = prec.prior, graph = adj.mat, constr = TRUE),
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

#' Note these match:

summary(inlafit)
summary(inlafitC)

#' However, if I change `constr = FALSE`, they no longer match even thought I think they
#' are fitting the same thing.

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

data$IDf <- factor(data$ID, 1:max(as.integer(data$ID)))

data2 <- filter(data, ID %in% 1:3, Year %in% 1973:1974) %>%
  droplevels()
adj.mat2 <- Matrix::spMatrix(3, 3, c(1, 2, 2, 3), c(2, 1, 3, 2), rep(1L, 4))
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

sdr <- TMB::sdreport(obj)
summary(sdr, "all")


#' Hyper parameters

c(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))
summary(sdr, "fixed")

#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr, "random")[1, ]

#' Random effects

inlafit$summary.random[[1]]

plot(inlafit$summary.random[[1]][ , 2], summary(sdr, "random")[-1, 1])
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][ , 3], summary(sdr, "random")[-1, 2])
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

data$IDf <- factor(data$ID, 1:max(as.integer(data$ID)))

data2 <- filter(data, ID %in% 1:3, Year %in% 1973:1974) %>%
  droplevels()
adj.mat2 <- Matrix::spMatrix(3, 3, c(1, 2, 2, 3), c(2, 1, 3, 2), rep(1L, 4))
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

tmbfit <- nlminb(obj$par, obj$fn, obj$gr,
                 control = list(iter.max = 1000,
                                eval.max = 1000))

tmbfit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

sdr <- TMB::sdreport(obj)
summary(sdr, "all")

prec.prior <- list(prec = list(param = c(0.001, 0.001)))


#' Hyper parameters

c(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))
summary(sdr, "fixed")

#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr, "random")[1, ]

#' Random effects

inlafit$summary.random[[1]]

plot(inlafit$summary.random[[1]][ , 2], summary(sdr, "report")[ , 1])
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][ , 3], summary(sdr, "report")[ , 2])
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

data$IDf <- factor(data$ID, 1:max(as.integer(data$ID)))
data$Yearf <- factor(data$Year, unique(data$Year))

data2 <- filter(data, ID %in% 1:3, Year %in% 1973:1974) %>%
  droplevels()
adj.mat2 <- Matrix::spMatrix(3, 3, c(1, 2, 2, 3), c(2, 1, 3, 2), rep(1L, 4))
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

tmbfit <- nlminb(obj$par, obj$fn, obj$gr,
                 control = list(iter.max = 1000,
                                eval.max = 1000))

sdr <- TMB::sdreport(obj)
summary(sdr, "all")

prec.prior <- list(prec = list(param = c(0.001, 0.001)))

inlafit <- inla(Observed ~ f(Year, model = "rw1", hyper = prec.prior),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))

summary(inlafit)

#' Hyper parameters

c(inlafit$misc$theta.mode,
  sqrt(diag(inlafit$misc$cov.intern)))

summary(sdr, "fixed")

#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr, "random")[1, ]

#' Random effects

inlafit$summary.random[[1]]

plot(inlafit$summary.random[[1]][ , 2], summary(sdr, "report")[ , 1])
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][ , 3], summary(sdr, "report")[ , 2])
abline(0, 1, col = "red")


#' #### INLA with custom Cmatrix

diagval <- INLA:::inla.set.f.default()$diagonal

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

data$IDf <- factor(data$ID, 1:max(as.integer(data$ID)))
data$Yearf <- factor(data$Year, unique(data$Year))

R_space <- diag(rowSums(adj.mat)) - adj.mat
R_space_adj <- R_space + Matrix::Diagonal(ncol(R_space), 1e-6)

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), 1e-6)

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

tmbfit <- nlminb(obj$par, obj$fn, obj$gr,
                 control = list(iter.max = 1000,
                                eval.max = 1000))

sdr <- TMB::sdreport(obj)
summary(sdr, "all")

prec.prior <- list(prec = list(param = c(0.001, 0.001)))

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

summary(sdr, "fixed")

#' Intercept

inlafit$summary.fixed[ , 1:2]
summary(sdr, "random")[1, ]

#' Random effects

sdrsum <- summary(sdr, "all")

plot(inlafit$summary.random[[1]][ , 2], sdrsum[rownames(sdrsum) == "u_time", 1],
     main = "f(Year): mean")
abline(0, 1, col = "red")

plot(inlafit$summary.random[[1]][ , 3], sdrsum[rownames(sdrsum) == "u_time", 2],
     main = "f(Year): sd")
abline(0, 1, col = "red")


plot(inlafit$summary.random[[2]][ , 2], sdrsum[rownames(sdrsum) == "u_space", 1],
     main = "f(area): mean")
abline(0, 1, col = "red")

plot(inlafit$summary.random[[2]][ , 3], sdrsum[rownames(sdrsum) == "u_space", 2],
     main = "f(area): sd")
abline(0, 1, col = "red")

#' #### INLA with custom Cmatrix

diagval <- INLA:::inla.set.f.default()$diagonal

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

data$IDf <- factor(data$ID, 1:max(as.integer(data$ID)))
data$Yearf <- factor(data$Year, unique(data$Year))

R_space <- diag(rowSums(adj.mat)) - adj.mat
R_space_adj <- R_space + Matrix::Diagonal(ncol(R_space), 1e-6)

R_time <- Matrix::sparseMatrix(1:19, 1:19, x = rep(1L, 19))

tmbdata <- list(y = data$Observed,
                E = data$Expected,
                Z_space_time = Matrix::sparse.model.matrix(~0 + IDf:Yearf, data),
                R_time = R_time,
                R_space = R_space)

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

sdr <- TMB::sdreport(obj)
summary(sdr, "all")

prec.prior <- list(prec = list(param = c(0.001, 0.001)))

inlafit <- inla(Observed ~
                  f(as.integer(IDf), model = "besag", hyper = prec.prior,
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

summary(sdr, "fixed")


#' ### Fixed effects (Intercept)

inlafit$summary.fixed[ , 1:2]
summary(sdr, "random")[1, , drop = FALSE]


#' ### Random effects mean and standard deviation

plot(inlafit$summary.random[[1]][,2], summary(sdr, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1)

plot(inlafit$summary.random[[1]][,3], summary(sdr, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1)


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
#' The weird thing about the results here is that the results here is that INLA returns
#' different results for the precision hyperparameter. But the estimates and standard
#' deviations for the fixed and random effects are the same.
#'
#' The results are not reproducible in TMB though with the same constraints.
#' 

inlafit <- inla(Observed ~ 
                  f(as.integer(IDf), model = "iid", hyper = prec.prior,
                    graph = adj.mat, diagonal = 1e-6, 
                    group = ID.Year, control.group = list(model = "rw1", scale.model = TRUE)),
                data = data, E = Expected, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute = list(config = TRUE))

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
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time,
                     hyper = prec.prior, diagonal = 1e-6),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))


#' These DON'T match:

summary(inlafit)
summary(inlafitC)

#' Hyper parameters

inlafit$internal.summary.hyperpar[, 1:2]
inlafitC$internal.summary.hyperpar[, 1:2]

#' But the random and fixed effect estimates both match...
#'
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

sdr <- TMB::sdreport(obj)
summary(sdr, "all")


#' ### Hyper parameter comparison
#'
#' The standard devation from TMB is somewhat larger than INLA.

summary(sdr, "fixed")

cbind("mean" = inlafit$misc$theta.mode,
      "se" = sqrt(diag(inlafit$misc$cov.intern)))

cbind("mean" = inlafitC$misc$theta.mode,
      "se" = sqrt(diag(inlafitC$misc$cov.intern)))



#' ### Fixed effects (Intercept)

summary(sdr, "random")[1, , drop = FALSE]
inlafit$summary.fixed[ , 1:2]
inlafitC$summary.fixed[ , 1:2]



#' ### Random effects mean and standard deviation

plot(inlafit$summary.random[[1]][,2], summary(sdr, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1)

plot(inlafit$summary.random[[1]][,3], summary(sdr, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1)






#' ## Model 7: ICAR space x RW1
#'
#' Change the time trend to a RW1. The TMB code remains the same. The only change is the
#' structure matrix for the time trend. If we were using a scaled parameterisation, we
#' would need to adjust the penalty for the additional rank deficiency of the structure
#' matrix.
#' 
#' These results show the same pattern as above where the two versions of the INLA
#' model give different hyper-parameter estimates, but the same fixed and random effects
#' 
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

#' #### INLA with custom Cmatrix

data$area.year <- interaction(data$IDf, data$Year)
data$id.area.year <- as.integer(data$area.year)

R_space <- diag(rowSums(adj.mat)) - adj.mat

D_time <- diff(diag(length(levels(data$Yearf))), differences = 1)
R_time <- Matrix::Matrix(t(D_time) %*% D_time)
R_time_adj <- R_time + Matrix::Diagonal(ncol(R_time), 1e-6)

R_space_time <- kronecker(R_time, R_space)

## Aconstr <- t(model.matrix(~0+factor(ID.Year), data[order(data$id.area.year), ]))
Aconstr <- t(apply(cbind(0:18, 1, 18:0)*32, 1, rep, x = c(0, 1, 0)))

inlafitC <- inla(Observed ~
                   f(id.area.year, model = "generic0", Cmatrix = R_space_time,
                     hyper = prec.prior, diagonal = 1e-6,
                     extraconstr = list(A = Aconstr, e = numeric(nrow(Aconstr)))),
                 data = data, E = Expected, family = "poisson",
                 control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                 control.compute = list(config = TRUE))


#' These DON'T match:

summary(inlafit)
summary(inlafitC)

#' Hyper parameters

inlafit$internal.summary.hyperpar[, 1:2]
inlafitC$internal.summary.hyperpar[, 1:2]

#' But the random and fixed effect estimates both match...
#'
#' Fixed effects

inlafit$summary.fixed
inlafitC$summary.fixed


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

sdr <- TMB::sdreport(obj)
summary(sdr, "all")


#' ### Hyper parameter comparison
#'
#' The standard devation from TMB is somewhat larger than INLA.

summary(sdr, "fixed")

cbind("mean" = inlafit$misc$theta.mode,
      "se" = sqrt(diag(inlafit$misc$cov.intern)))

cbind("mean" = inlafitC$misc$theta.mode,
      "se" = sqrt(diag(inlafitC$misc$cov.intern)))



#' ### Fixed effects (Intercept)

summary(sdr, "random")[1, , drop = FALSE]
inlafit$summary.fixed[ , 1:2]
inlafitC$summary.fixed[ , 1:2]



#' ### Random effects mean and standard deviation

plot(inlafit$summary.random[[1]][,2], summary(sdr, "report")[,1],
     xlab = "INLA", ylab = "TMB", main = "Random effect point estimates")
abline(0, 1)

plot(inlafit$summary.random[[1]][,3], summary(sdr, "report")[,2],
     xlab = "INLA", ylab = "TMB", main = "Random effect standard deviation")
abline(0, 1)





