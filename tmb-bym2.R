#' ---
#' title: INLA sparsity preserving BYM2 implementation in TMB
#' author: Jeff Eaton (jeffrey.eaton@imperial.ac.uk)
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 2
#' ---
#'
#' # INLA parameterisation for the BYM2 model
#' 
#' We frequently the _BYM2_ model for spatially structured random effects. The 
#' BYM2 model is a reparameterisation of the _Besag-York-Mollie (BYM)_ model to
#' decompose a spatial random effect $x$ into a spatially structured component
#' $u$ and a an IID component $v$ defined by a marginal standard deviation
#' $\sigma > 0$ and proportion $\phi \in [0,1]$ specifying the weight to
#' each component:
#'
#' $$ \mathbf{x} = \sigma \cdot (\sqrt{\phi}\cdot \mathbf{u} + \sqrt{1-\phi}\cdot \mathbf{v}) $$
#' $$ \mathbf{u} \sim N(0, \mathbf{Q}^{-1}) $$
#' $$ \mathbf{v} \sim N(0, \mathbf{I}) $$
#'
#' where the matrix $\mathbf{Q}$ is the _scaled_ ICAR structure matrix such that
#' $\mathrm{Var}(u_i) \approx 1$.
#'
#' In most TMB implementations, we parameterise the BYM2 model in terms of the
#' component random effects $\mathbf{u}$ and $\mathbf{v}$, as illustrated in the
#' [tutorial on intrinsic area models in Stan](https://mc-stan.org/users/documentation/case-studies/icar_stan.html#bym2-improving-the-parameterization-of-the-besag-york-and-mollie-model)
#' by Mitzi Morris.
#'
#' The INLA implementation of the BYM2 model directly parameterises the spatial
#' random effect $\mathbf{x}$, which has covariance matrix
#' $$\mathrm{Var}(\mathbf{x} | \sigma, \phi) = \sigma^2\cdot\left(\phi\cdot\mathbf{Q}^{-1} + (1-\phi)\cdot\mathbf{I}\right).$$
#'
#' The inverse of $\mathrm{Var}(\mathbf{x} | \sigma, \phi)$, the precision matrix,
#' is no longer sparse. To retain sparsity of the precision matrix, INLA parameterises
#' the joint distribution of $(\mathbf{x}, \mathbf{u})$, described in Section 3.4 of
#' [Reibler _et al._](https://arxiv.org/pdf/1601.01180.pdf). 
#' Conditioning $\pi(\mathbf{x}, \mathbf{u}) = \pi(\mathbf{x} | \mathbf{u}) \pi(\mathbf{u})$
#' is used to derive that $\pi(\mathbf{x}, \mathbf{u})$ is a multivariate normal distribution
#' with mean $\mathbf{0}$ and precision matrix
#'
#' $$\left(\begin{array}{cc}
#' \frac{1}{\sigma^2 (1-\phi)} \mathbf{I} & -\frac{\sqrt{\phi}}{\sigma (1-\phi)} \mathbf{I} \\
#' -\frac{\sqrt{\phi}}{\sigma (1-\phi)} \mathbf{I} & \mathbf{Q} + \frac{\phi}{1-\phi} \mathbf{I}
#' \end{array}\right).$$
#'
#' The normalising constant is apparent from the product of the normalising constants
#' for the terms $\pi(\mathbf{x} | \mathbf{u})$ and $\pi(\mathbf{u})$:
#' 
#' $$ (2\pi)^{-(n + \mathrm{rank}(\mathbf{Q}))/2} \cdot \left(\sigma \sqrt{1-\phi}\right)^{-n} \cdot |\mathbf{Q}|^{1/2}. $$
#'
#' The terms $(2\pi)^{-(n + \mathrm{rank}(\mathbf{Q}))/2}$ and $|\mathbf{Q}|^{1/2}$
#' typically may be omitted.
#' Further details of the derivation are on page 3 of the [notes accompanying the
#' INLA documentation](https://github.com/hrue/r-inla/blob/devel/r-inla.org/doc/latent/bym2-details.pdf).
#'
#' ## Benefit of directly parameterising $\mathbf{x}$
#'
#' TMB, like INLA, implements a two step optimisation process. The 'outer' optimisation
#' step maximises the marginal posterior of the hyperparameters (referred to as 'fixed'
#' parmaeters in TMB). For each 'outer' optimisation step, an 'inner' optimisation
#' maximises the latent field parameters ('random' parameters in TMB parlance) 
#' conditional on the hyperparameters, which is used for Laplace approximation to the
#' marginal posterior of the hyperparameters.
#' 
#' The starting values for each inner optimisation are the final optimised values from
#' the previous iteration. Using the conventional parameterisation from the Morris
#' tutorial, the linear predictor for an observation depends on the values of the
#' hyperparameters, for example
#' $$\mu_i = \beta_0 + \sigma\cdot\left(\sqrt{\phi}\cdot u_i + \sqrt{1-\phi}\cdot v_i\right).$$
#' Thus a large step for the hyperparameters will affect the the value for the linear
#' predictor $\mu_i$, potentially moving it far away from the optimal parameters.
#'
#' If the random effect $x_i$ is directly parameterised as in the INLA implementation,
#' for example $\mu_i = \beta_0 + x_i,$ the starting value for $\mu_i$ is unchanged
#' by the updated hyperparamters. This change in parameterisation does not change
#' the model (_i.e._ the inference for hyperparameters does not change), but
#' the inner optimisation step may be more efficient by starting at parameter values
#' closer to the final values.
#'
#' ## Implementing the sparsity preserving parameterisation in TMB
#'
#' INLA implements the BYM2 model parameterisation as a single vector 
#' $\mathbf{w} = (\mathbf{x}, \mathbf{u})$ of length $2n$, and constructs the
#' full $2n \times 2n$ sparse precision matrix described above. Coding the
#' sparse precision matrix and parsing the vector $\mathbf{w}$ in the TMB model
#' template is an unnecessary hassle.
#'
#' Instead, it is more convenient to decompose the block matrix multiplication:
#'
#' \begin{align*}
#' \log \pi(\mathbf{x}, \mathbf{u} | \sigma,\phi) &= -\frac{n+\mathrm{rank}(\mathbf{Q})}{2}
#'    \log(2\pi) - n\log(\sigma) - \frac{n}{2} \log(1-\phi) + \frac{1}{2}\log|\mathbf{Q}| - 
#'   \frac{1}{2} \left(\begin{array}{cc}\mathbf{x}' & \mathbf{u}'\end{array}\right)
#'      \left(\begin{array}{cc}
#' \frac{1}{\sigma^2 (1-\phi)} \mathbf{I} & -\frac{\sqrt{\phi}}{\sigma (1-\phi)} \mathbf{I} \\
#' -\frac{\sqrt{\phi}}{\sigma (1-\phi)} \mathbf{I} & \mathbf{Q} + \frac{\phi}{1-\phi} \mathbf{I}
#' \end{array}\right)
#'      \left(\begin{array}{c} \mathbf{x} \\ \mathbf{u}\end{array}\right) \\
#'   &= C - \frac{1}{2}\left[\mathbf{x}' \left(\frac{1}{\sigma^2 (1-\phi)} \mathbf{I}\right)\mathbf{x} -
#'      2\mathbf{x}'\left(\frac{\sqrt{\phi}}{\sigma (1-\phi)} \mathbf{I}\right)\mathbf{u} +
#'      \mathbf{u}'\left(\mathbf{Q} + \frac{\phi}{1-\phi} \mathbf{I}\right) \mathbf{u}\right] \\
#'   &= -\frac{n+\mathrm{rank}(\mathbf{Q})}{2}
#'    \log(2\pi) - n\log(\sigma) - \frac{n}{2} \log(1-\phi) + \frac{1}{2}\log|\mathbf{Q}| -
#'    \frac{\mathbf{x}'\mathbf{x}}{2\sigma^2(1-\phi)} +
#'      \frac{\sqrt{\phi}\mathbf{x}'\mathbf{u}}{\sigma (1-\phi)} -
#'      \frac{\phi \mathbf{u}'\mathbf{u}}{2\cdot(1-\phi)} -
#'      \frac{1}{2}\mathbf{u}'\mathbf{Q}\mathbf{u}
#' \end{align*}
#'
#' The constant terms $-\frac{n+\mathrm{rank}(\mathbf{Q})}{2}\log(2\pi)$ and $\frac{1}{2}\log|\mathbf{Q}|$
#' typically need not be computed.
#' 
#' # Example: Scottish Lip Cancer dataset
#'
#' ## Preliminaries
#'
#' Load libraries and utility functions.

##+ message = FALSE, results = "hide"
library(tidyverse)
library(INLA)
library(TMB)

#' The function `tmb_compile_and_load()` is a utility function that accepts a TMB
#' model as a string `code`, compiles and loads the model, and return a path to
#' the DLL.

tmb_compile_and_load <- function(code) {
  f <- tempfile(fileext = ".cpp")
  writeLines(code, f)
  TMB::compile(f)
  dyn.load(TMB::dynlib(tools::file_path_sans_ext(f)))
  basename(tools::file_path_sans_ext(f))
}


#' ## Scottish Lip Cancer Data
#' 
#' Use the version of Scottish lip cancer dataset used by Mitzi Morris for the
#' [tutorial on ICAR models in Stan](https://mc-stan.org/users/documentation/case-studies/icar_stan.html).

source("https://raw.githubusercontent.com/stan-dev/example-models/master/knitr/car-iar-poisson/scotland_data.R")
scotlip <- data

#' The list `scotlip` consists of items:
#'
#' * `N`: number of regions,
#' * `y`: observed counts of lip cancer cases per county,
#' * `E`: the expected number of cases, used as an offset,
#' * `x`: continuous covariate representing the proportion of the population employed in agriculture, fishing, or forestry (AFF),
#' * `adj`: a vector of region ids for adjacent regions,
#' * `weights`: weight for each edge (all `1`s),
#' * `num`: the number of neighbors for each region, used to split `adj`.
#'
#' Scale the covariate `x` as in Morris tutorial.

scotlip$x_scaled <- 0.1 * scotlip$x

#' Parse the adjacency list into an adjacency matrix.

nblist <- cbind(rep(seq_len(scotlip$N), times = scotlip$num),
                scotlip$adj)
adj <- matrix(0, nrow = scotlip$N, ncol = scotlip$N)
adj[nblist] <- 1


#' Construct the scaled structure matrix for ICAR model precision.

Q <- diag(rowSums(adj)) - adj
Q_scaled <- inla.scale.model(Q, constr = list(A = matrix(1, ncol = ncol(Q)), e = 0))


#' Prepare a list of TMB model inputs.

tmbdata <- list(region = 1:scotlip$N,
                y = scotlip$y,
                x = scotlip$x_scaled,
                E = scotlip$E,
                Q = Q_scaled)


#' ## 1. Morris parameterisation: $\mathbf{u}$ and $\mathbf{v}$
#'
#' The TMB model template `mod1` below implements the model
#' $$ y \sim \mathrm{Poisson}(\eta_i) $$
#' $$ \log \eta_i = \beta_0 + \beta_1 \cdot x + b_i + \log E_i $$
#' $$ \mathbf{b} = \sigma \cdot \left(\sqrt{\phi}\cdot \mathbf{u} + \sqrt{1-\phi}\cdot \mathbf{v}\right) $$
#' $$ \mathbf{u} \sim N(0, \mathbf{Q}^{-1}) $$
#' $$ \mathbf{v} \sim N(0, \mathbf{I}) $$
#' $$ \sigma \sim N^{+}(0, 1) $$
#' $$ \phi \sim \mathrm{Beta}(0.5, 0.5) $$
#' 

mod1 <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR area model

  Type val(0);

  PARAMETER(beta0);  // intercept
  PARAMETER(beta1);  // slope
  // beta0 ~ 1
  // beta1 ~ 1

  PARAMETER(log_sigma); // marginal standard deviation
  Type sigma(exp(log_sigma));
  val -= dnorm(sigma, Type(0.0), Type(1.0), true) + log_sigma;

  PARAMETER(logit_phi);
  Type phi(invlogit(logit_phi));
  val -= log(phi) +  log(1 - phi);  // change of variables: logit_phi -> phi
  val -= dbeta(phi, Type(0.5), Type(0.5), true);

  PARAMETER_VECTOR(u); // spatially correlated component
  val -= Type(-0.5) * (u * (Q * u)).sum();

  // soft sum-to-zero constraint
  val -= dnorm(u.sum(), Type(0.0), Type(0.001) * u.size(), true); 

  PARAMETER_VECTOR(v); // unstructured component
  val -= dnorm(v, Type(0.0), Type(1.0), true).sum();

  // combined spatial effect
  vector<Type> b(sigma * (sqrt(phi) * u + sqrt(1 - phi) * v));

  vector<Type> mu(beta0 + beta1 * x + b + log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(sigma);
  ADREPORT(phi);
  ADREPORT(b);
  ADREPORT(mu);

  return val;
}
'

#' Compile and fit the TMB model.
#' 
##+ message = FALSE, results = "hide"
dll1 <- tmb_compile_and_load(mod1)

#' Initial values for parameters.

tmbpar1 <- list(beta0 = 0,
                beta1 = 0,
                log_sigma = 0,
                logit_phi = 0,
                u = numeric(scotlip$N),
                v = numeric(scotlip$N))

#' Create TMB object and optimise parameters.

##+ message = FALSE, results = "hide"
obj1 <- TMB::MakeADFun(data = tmbdata,
                       parameters = tmbpar1,
                       random = c("beta0", "beta1", "u", "v"),
                       DLL = dll1)

tmbfit1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

sdr1 <- TMB::sdreport(obj1)

#' Estimates for intercept and slope and hyperparameters.

summary(sdr1, "all") %>%
  .[rownames(.) %in% c("log_sigma", "logit_phi", "beta0", "beta1"), ]

#' ## 2. INLA parameterisation: $\mathbf{b}$ and $\mathbf{u}$
#' 
#' The TMB model template `mod2` below implements the INLA parameterisation
#' of the same model described in Section 3.4 of Reibler _et al._:
#' $$ y \sim \mathrm{Poisson}(\eta_i) $$
#' $$ \log \eta_i = \beta_0 + \beta_1 \cdot x + b_i + \log E_i $$
#' $$ (\mathbf{b}, \mathbf{u}) \sim N(0, \mathbf{\Sigma}(\sigma, \phi)) $$
#' $$ \sigma \sim N^{+}(0, 1) $$
#' $$ \phi \sim \mathrm{Beta}(0.5, 0.5) $$
#'
#' where the precision matrix $\mathbf{\Sigma}^{-1}(\sigma, \phi)$ is as defined above.
#' 

mod2 <- '
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_VECTOR(E);

  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR area model

  Type val(0);

  PARAMETER(beta0);  // intercept
  PARAMETER(beta1);  // intercept
  // beta0 ~ 1
  // beta1 ~ 1

  PARAMETER(log_sigma); // marginal standard deviation
  Type sigma(exp(log_sigma));
  val -= dnorm(sigma, Type(0.0), Type(1.0), true) + log_sigma;

  PARAMETER(logit_phi);
  Type phi(invlogit(logit_phi));
  val -= log(phi) +  log(1 - phi);  // change of variables: logit_phi -> phi
  val -= dbeta(phi, Type(0.5), Type(0.5), true);

  PARAMETER_VECTOR(b); // combined spatial effect

  PARAMETER_VECTOR(u); // spatially correlated component
  val -= dnorm(u.sum(), Type(0.0), Type(0.001) * u.size(), true); // soft sum-to-zero constraint

  // -(n + rank(Q)) / 2 * log(2*pi) and log|Q| / 2 terms omitted
  val -= -0.5 * b.size() * (2 * log_sigma + log(1 - phi));          // normalising constant
  val -= -0.5 * (b * b).sum() / (sigma * sigma * (1 - phi));
  val -= (b * u).sum() * sqrt(phi) / (sigma * (1 - phi));
  val -= -0.5 * (u * u).sum() * phi / (1 - phi);
  val -= -0.5 * (u * (Q * u)).sum();

  vector<Type> mu(beta0 + beta1 * x + b + log(E));
  val -= dpois(y, exp(mu), true).sum();

  ADREPORT(sigma);
  ADREPORT(phi);
  ADREPORT(mu);

  return val;
}
'

#' Compile and fit model 2.

##+ message = FALSE, results = "hide"
dll2 <- tmb_compile_and_load(mod2)

tmbpar2 <- list(beta0 = 0,
                beta1 = 0,
                log_sigma = 0,
                logit_phi = 0,
                b = numeric(scotlip$N),
                u = numeric(scotlip$N))

obj2 <- TMB::MakeADFun(data = tmbdata,
                       parameters = tmbpar2,
                       random = c("beta0", "beta1", "b", "u"),
                       DLL = dll2)

tmbfit2 <- nlminb(obj2$par, obj2$fn, obj2$gr)

sdr2 <- TMB::sdreport(obj2)

#' ## Compare estimates
#' 
#' Estimates for intercept and slope and hyperparameters match exactly with the
#' Morris parameterisation.

summary(sdr1, "all") %>%
  .[rownames(.) %in% c("log_sigma", "logit_phi", "beta0", "beta1"), ]

summary(sdr2, "all") %>%
  .[rownames(.) %in% c("log_sigma", "logit_phi", "beta0", "beta1"), ]

#' The scatter plot of point estimates and posterior standard deviation for the
#' random effects $b_i$ and $u_i$ show these are also identical.
#' 
##+ echo = FALSE, fig.align = "center", fig.height = 6, fig.width = 6
est <- summary(sdr1, "all") %>%
  data.frame(par = rownames(.), ., version = "Version 1 (Morris)") %>%
  group_by(par) %>%
  mutate(region = row_number()) %>%
  bind_rows(
    summary(sdr2, "all") %>%
    data.frame(par = rownames(.), ., version = "Version 2 (INLA)") %>%
    group_by(par) %>%
    mutate(region = row_number()) 
  )

est %>%
  pivot_longer(c(Estimate, Std..Error)) %>%
  pivot_wider(names_from = version) %>%
  filter(par %in% c("b", "u")) %>%
  ggplot(aes(`Version 1 (Morris)`, `Version 2 (INLA)`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_point() +
  facet_wrap(~ par + name, nrow = 2, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())

#' ## Check TMB results with INLA
#'
#' The TMB models above conduct the same inference as INLA using the options `control.inla = list(strategy = "gaussian", int.strategy = "eb")`. 
#' 
#' 

##+ warning = FALSE, message = FALSE
priors <- list(phi = list(prior = "logitbeta", params = c(0.5, 0.5)),
               prec = list(prior = "logtnormal", params = c(0, 1.0)))

inla_formula <- y ~ 1 + x +
  f(region, model = "bym2", graph = adj, hyper = priors, constr = TRUE)

inla_eb <- inla(inla_formula, family = "poisson", E = E, 
                data = tmbdata[c("region", "y", "x", "E")],
                control.fixed = list(prec = 1/25, prec.intercept = 1/25),
                control.inla = list(strategy = "gaussian", int.strategy = "eb"))

#' Estimates and standard errors for the fixed effects and hyperparameters from TMB
#' are close to the mode and sd reported by INLA. The log precision is -2 times the 
#' `log_sigma` parameter in the TMB models.

inla_eb$internal.summary.hyperpar
inla_eb$summary.fixed

#' Point estimates for random effect estimates are exactly the same
#' for TMB and the INLA empirical Bayes estimates. The posterior
#' standard deviations are highly correlated but slightly larger for
#' the TMB estimates than INLA.

##+ echo = FALSE, fig.align = "center", fig.height = 6, fig.width = 6
inla_eb$summary.random[[1]] %>%
  mutate(par = rep(c("b", "u"), each = nrow(.) / 2),
         region = rep(seq_len(nrow(.)/2), times = 2),
         version = "INLA empirical Bayes") %>%
  rename(Estimate = mode, Std..Error = sd) %>%
  select(all_of(names(est))) %>%
  bind_rows(
    filter(est, version == "Version 2 (INLA)")
  ) %>%
  pivot_longer(c(Estimate, Std..Error)) %>%
  pivot_wider(names_from = version) %>%
  filter(par %in% c("b", "u")) %>%
  ggplot(aes(`Version 2 (INLA)`, `INLA empirical Bayes`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_point() +
  facet_wrap(~ par + name, nrow = 2, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())

#' # Full Bayesian inference
#'
#' The point estimates from the TMB and INLA models above are slightly different
#' from (those reported in the Morris tutorial)[https://mc-stan.org/users/documentation/case-studies/icar_stan.html#bym2-improving-the-parameterization-of-the-besag-york-and-mollie-model]
#' using full Bayesian inference in Stan or INLA.
#' 
#' These are reproduced by changing the strategy to `control.inla = list(strategy="laplace")`.

inla_full <- inla(inla_formula, family = "poisson", E = E, 
                  data = tmbdata[c("region", "y", "x", "E")],
                  control.fixed = list(prec = 1/25, prec.intercept = 1/25),
                  control.inla = list(strategy="laplace"))

inla_full$summary.fixed
inla_full$internal.summary.hyperpar

#' The package [`tmbstan`](https://github.com/kaskr/tmbstan) can be used to sample from
#' the posterior distribution with the TMB model objects created above. 
#' 

##+ warning = FALSE, message = FALSE
library(tmbstan)

stanfit1 <- tmbstan(obj1, refresh = 0)
stanfit2 <- tmbstan(obj2, refresh = 0)

print(stanfit1, digits = 3, par = c("beta0", "beta1", "log_sigma", "logit_phi", "lp__"))
print(stanfit2, digits = 3, par = c("beta0", "beta1", "log_sigma", "logit_phi", "lp__"))

#' The estimates from `stanfit1`, using the Morris tutorial parameterisation, match closely
#' to the full INLA results.
#'
#' For the second parameterisation, the Rhat estimate for the `logit_phi` parameter
#' is large and the number of effective samples is low, indicating poor convervenge
#' for this parameter. This highlights that direct parameterisation of $\mathbf{b}$
#' is not optimal for HMC inference, at least in this case.
#' 
