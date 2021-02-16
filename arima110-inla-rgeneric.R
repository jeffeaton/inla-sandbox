library(INLA)
library(Matrix)
library(ggplot2)
library(dplyr)
library(forcats)

#' ## INLA AR(1) example using rgeneric.
#' From https://inla.r-inla-download.org/r-inla.org/doc/vignettes/rgeneric.pdf.
#'

inla.rgeneric.ar1.model <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                            "log.prior", "quit"),
                                    theta = NULL) {

  ## for reference and potential storage for objects to
  ## cache, this is the environment of this function
  ## which holds arguments passed as `...` in
  ## `inla.rgeneric.define()`.

  interpret.theta = function() {
    return(list(prec = exp(theta[1L]),
                rho = 2 * exp(theta[2L])/(1 + exp(theta[2L])) - 1))
  }
  envir = parent.env(environment())
  
  graph = function(){ return (Q()) }
  
  Q = function() {
    p = interpret.theta()
    i = c(1L, n, 2L:(n - 1L), 1L:(n - 1L))
    j = c(1L, n, 2L:(n - 1L), 2L:n)
    x = p$prec/(1 - p$rho^2) *
      c(1L, 1L, rep(1 + p$rho^2, n - 2L),
        rep(-p$rho, n - 1L))
    return (sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE))
  }

  mu = function() { return(numeric(0)) }

  log.norm.const = function() {
    p = interpret.theta()
    prec.i = p$prec / (1.0 - p$rho^2)
    val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.i)) +
      0.5 * log(1.0 - p$rho^2)
    return (val)
  }
  
  log.prior = log.prior = function() {
    p = interpret.theta()
    val = dgamma(p$prec, shape = 1, rate = 1, log=TRUE) + theta[1L] +
      dnorm(theta[2L], mean = 0, sd = 1, log=TRUE)
    return (val)
  }

  initial = function() {
    return (rep(1, 2))
  }

  quit = function() { return (invisible()) }

  
  ## sometimes this is useful, as argument 'graph' and 'quit'
  ## will pass theta=numeric(0) (or NULL in R-3.6...) as
  ## the values of theta are NOT
  ## required for defining the graph. however, this statement
  ## will ensure that theta is always defined.
  if (!length(theta)) theta = initial()

  val = do.call(match.arg(cmd), args = list())
  return (val)
}


n <- 100
rho <- 0.9
x <- arima.sim(n, model = list(ar = rho)) * sqrt(1-rho^2)
y <- x + rnorm(n, sd = 0.1)
model <- inla.rgeneric.define(inla.rgeneric.ar1.model, n = n)
formula <- y ~ -1 + f(idx, model = model)
r <- inla(formula, data = data.frame(y, idx = 1:n))

summary(r)


fformula = y ~ -1 +
  f(idx, model = "ar1",
    hyper = list(prec = list(prior = "loggamma", param = c(1,1)),
                 rho = list(prior = "normal", param = c(0,1))))

rr = inla(fformula, data = data.frame(y, idx = 1:n))


plot(inla.smarginal(rr$marginals.hyperpar[[2]]),
     type="l", lwd=5, col="red", xlab="stdev", ylab="density")
lines(inla.tmarginal(exp, r$internal.marginals.hyperpar[[2]]),
      col="yellow")


plot(inla.smarginal(rr$marginals.hyperpar[[3]]),
     type="l", lwd=5, col="red", xlab="rho", ylab="density")
lines(inla.tmarginal(function(x) 2*exp(x)/(1+exp(x))-1,
                     r$internal.marginals.hyperpar[[3]]),
      col="yellow")

round(rbind(native = rr$cpu.used,
            rgeneric = r$cpu.used), digits = 3)


#' ## ARIMA(1, 1, 0)
#'
#' $\alpha ~ ARIMA(1, 1, 0)$ is equivalent to $\gamma ~ AR(1)$
#' where $\gamma = D\alpha$ and D is the first order difference
#' matrix.
#'
#' Let Q~ be the AR(1) precision matrix. Then
#' $\gamma' Q~ \gamma = \alpha' D' Q~ D \alpha$.
#'
#' Thus precision matrix for ARIMA(1, 1, 0) is simly $Q = D' Q~ D$.
#'
#' det(Q) = n * det(Q~).  I don't know why yet...
#'
#' The normalising constant has coefficient n-1.
#' 


inla.rgeneric.arima110.model <- function(cmd = c("graph", "Q", "mu", "initial",
                                                 "log.norm.const", "log.prior", "quit"),
                                         theta = NULL) {

  ## for reference and potential storage for objects to
  ## cache, this is the environment of this function
  ## which holds arguments passed as `...` in
  ## `inla.rgeneric.define()`.

  interpret.theta = function() {
    return(list(prec = exp(theta[1L]),
                rho = 2 * exp(theta[2L])/(1 + exp(theta[2L])) - 1))
  }
  envir = parent.env(environment())
  
  graph = function(){ return (Q()) }
  
  Q = function() {

    n.ar1 <- n-1  # dimension of AR1 on first order differences
    D1 = Matrix::diff(Matrix::Diagonal(n))  # first order difference operator
    
    p = interpret.theta()

    ## Define both uperr and lower tri because D'*Q*D calculation
    i = c(1L, n.ar1, 2L:(n.ar1 - 1L), 1L:(n.ar1 - 1L), 2L:n.ar1)
    j = c(1L, n.ar1, 2L:(n.ar1 - 1L), 2L:n.ar1, 1L:(n.ar1 - 1L))
    x = p$prec/(1 - p$rho^2) *
      c(1L, 1L, rep(1 + p$rho^2, n.ar1 - 2L),
        rep(-p$rho, n.ar1 - 1L),
        rep(-p$rho, n.ar1 - 1L))

    Q.ar1 <- sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE)
    
    return (t(D1) %*% Q.ar1 %*% D1)
  }

  mu = function() { return(numeric(0)) }

  log.norm.const = function() {
    p = interpret.theta()
    prec.i = p$prec / (1.0 - p$rho^2)
    val = (n - 1) * (- 0.5 * log(2*pi) + 0.5 * log(prec.i)) +
      0.5 * log(1.0 - p$rho^2) +
      log(n)  # + log(n) because det(Q) = n * det(Q.ar1)
    
    return (val)
  }
  
  log.prior = log.prior = function() {
    p = interpret.theta()
    val = dgamma(p$prec, shape = 1, rate = 1, log=TRUE) + theta[1L] +
      dnorm(theta[2L], mean = 0, sd = 1, log=TRUE)
    return (val)
  }

  initial = function() {
    return (rep(1, 2))
  }

  quit = function() { return (invisible()) }

  
  ## sometimes this is useful, as argument 'graph' and 'quit'
  ## will pass theta=numeric(0) (or NULL in R-3.6...) as
  ## the values of theta are NOT
  ## required for defining the graph. however, this statement
  ## will ensure that theta is always defined.
  if (!length(theta)) theta = initial()

  val = do.call(match.arg(cmd), args = list())
  return (val)
}


#' Simulate data

set.seed(1)

n <- 100
rho <- 0.8
c <- 0.3

#' The parameter `mean = <>` to `arima.sim()` passes a mean value to `rnorm()` used to
#' simulate the innovations. To obtain a process with marginal slope `c`, pass the
#' value `mean = c * (1 - rho) / sqrt(1-rho^2)`.

x <- arima.sim(n, model = list(order = c(1, 1, 0), ar = rho), mean = c * (1 - rho) / sqrt(1-rho^2))
x <- x[-1] * sqrt(1-rho^2)

y <- x + rnorm(n, sd = 0.1)

#' Define INLA ARIMA(1,1,0) model object
arima110_model <- inla.rgeneric.define(inla.rgeneric.arima110.model, n = n)
                                       
#' Fit INLA model to full data
formula <- y ~ idx_lin + f(idx, model = arima110_model, constr = TRUE)
fit <- inla(formula, data = data.frame(y, idx_lin = 1:n, idx = 1:n),
            control.predictor = list(compute = TRUE))

summary(fit)

fit$summary.fitted.values %>%
  mutate(idx = 1:n) %>%
  ggplot(aes(idx, mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  geom_point(aes(y = y), color = "red3")


#' Truncate last 25% projection data

y_trunc <- y
y_trunc[round(n*0.75):n] <- NA

fit_trunc <- inla(y_trunc ~ idx_lin + f(idx, model = arima110_model, constr = TRUE),
                  data = data.frame(y_trunc, idx_lin = 1:n, idx = 1:n),
                  control.predictor = list(compute = TRUE))

summary(fit_trunc)

fit_trunc$summary.fitted.values %>%
  mutate(idx = 1:n) %>%
  ggplot(aes(idx, mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  geom_point(aes(y = y), color = "red3")


#' ## Compare ARIMA(1,1,0) + c with other models
#'
#' * ARIMA(1, 1, 0) [no constant]
#' * RW(2)
#' * AR(2)
#' * AR(2) + c

fit_arima110 <- inla(y_trunc ~ f(idx,  model = arima110_model, constr = TRUE),
                     data = data.frame(y_trunc, idx_lin = 1:n, idx = 1:n),
                     control.predictor = list(compute = TRUE),
                     control.compute = list(config = TRUE))

fit_rw2 <- inla(y_trunc ~ f(idx, model = "rw2"),
                data = data.frame(y_trunc, idx_lin = 1:n, idx = 1:n),
                control.predictor = list(compute = TRUE))

fit_ar2 <- inla(y_trunc ~ f(idx, model = "ar", order = 2),
                data = data.frame(y_trunc, idx_lin = 1:n, idx = 1:n),
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE))

fit_ar2_c <- inla(y ~ idx_lin + f(idx, model = "ar", order = 2),
                  data = data.frame(y_trunc, idx_lin = 1:n, idx = 1:n),
                  control.predictor = list(compute = TRUE))

fitted_all <- bind_rows(
  mutate(fit_trunc$summary.fitted.values, idx = row_number(), model = "ARIMA(1,1,0) + c"),
  mutate(fit_arima110$summary.fitted.values, idx = row_number(), model = "ARIMA(1,1,0)"),
  mutate(fit_rw2$summary.fitted.values, idx = row_number(), model = "RW2"),
  mutate(fit_ar2$summary.fitted.values, idx = row_number(), model = "AR(2)"),
  mutate(fit_ar2_c$summary.fitted.values, idx = row_number(), model = "AR(2) + c")
) %>%
  mutate(model = fct_inorder(model))

fitted_all %>%
  ggplot(aes(idx, mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_ribbon(aes(fill = model, color = model), alpha = 0.2, linetype = "dotted", size = 0.25) +  
  geom_line(aes(color = model)) +
  geom_point(aes(x = idx, y = y), color = "red3", inherit.aes = FALSE,
             data = data.frame(idx = seq_along(y), y)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")
