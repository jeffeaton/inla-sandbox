#' ---
#' title: INLA internals
#' author: Jeff Eaton (jeffrey.eaton@imperial.ac.uk)
#' output:
#'   md_document:
#'     variant: markdown_github
#'   
#' ---
#' 
#' This document is to record the details of some of the internal specifications 
#' of INLA. The purpose of this is largely to document what INLA is doing for 
#' comparing model implementations in other software.
#' 
#' **Make sure to be using a version of INLA more recent than INLA_20.06.15 _testing_
#' version when the internal implementation of the RW1 and RW2 models was updated to
#' be consistent with manual model scaling.**
#' 

library(INLA)

#' # Data
#' 
#' Create some very simple space-time data for testing INLA. Data are for four
#' areas in which area 2 is connected to all other areas, areas 1 and 3 are 
#' connected, and area 4 is connected only to area 2. Three time points are 
#' simulated.

##+ data
set.seed(1)
data <- expand.grid(area = 1:4,
                    time = 1:3)
                    
data$y <- rpois(nrow(data), 2.5)

adj <- rbind(c(0, 1, 1, 0),
             c(0, 0, 1, 1),
             0,
             0)
adj <- adj + t(adj)
rownames(adj) <- colnames(adj) <- letters[1:4]

adj

#' Structure matrices for ICAR area effect and RW1 time effects.

##+ structure_mat

R_area <- diag(rowSums(adj)) - adj
R_time <- t(diff(diag(3))) %*% diff(diag(3))

R_area
R_time

#' # Details of INLA internals
#' 
#' This section documents details of how certain arguments and options are 
#' implemented by INLA such as how the small constant is added to the diagonal,
#' model scaling, and kronecker product for the `f(..., group = <>)` option.
#' 
#' The strategy for checking is to create a 'null' data set with no observations
#' and 'fit' the model with fixed values for the hyper parameters to recover
#' the Q matrix constructed by INLA.

datanull <- data[data$area == 1, ]
datanull$y <- NA

hyper_area <- list(prec = list(initial = log(1), fixed = TRUE))
hyper_time <- list(prec = list(initial = log(1), fixed = TRUE))

#' ## ICAR model 
#' 

fit1 <- inla(y ~ 0 +
               f(area, model = "besag", graph = adj, hyper = hyper_area),
             data = datanull, family = "poisson",
             control.inla = list(strategy = "gaussian", int.strategy = "eb"),
             control.fixed = list(mean.intercept = 0, prec.intercept = 1),
             control.compute = list(config = TRUE))

#' ### Extra constant added to diagonal
#' 
#' * The default value for the diagonal extra constant is `r INLA:::inla.set.f.default()$diagonal`.
#'   * This is ascertained from `INLA:::inla.set.f.default()$diagonal`.
#' * This is added *after* scaling the matrix by the precision parameter.
#'
#' To see this, in `fit1` the fixed value for the precision is 1.0 and in
#' `fit2` the value for the precision is 2.0. The added value along the
#' diagonal in both cases is 1e-5.
#' 

hyper_area2 <- list(prec = list(initial = log(2), fixed = TRUE))

fit2 <- inla(y ~ 0 +
               f(area, model = "besag", graph = adj, hyper = hyper_area2),
             data = datanull, family = "poisson",
             control.inla = list(strategy = "gaussian", int.strategy = "eb"),
             control.fixed = list(mean.intercept = 0, prec.intercept = 1),
             control.compute = list(config = TRUE))

fit1$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]
fit2$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]

#' ### Application of `f(..., scale.model = TRUE)`
#' 
#' When argument `scale.model = TRUE`, the precision matrix is scaled so that
#' the generalised variance is 1. 
#' 
#' * A sum-to-zero constraint is applied in the model scaling. (Likely a 
#'   different constraint is applied if there are multiple connected components).
#' * No constant is added to the diagonal before model scaling.
#' 

fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                scale.model = TRUE),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

fit$misc$configs$constr

fit$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]

inla.scale.model(R_area, constr = list(A = matrix(1, ncol = 4), e = 0))

#'
#' * Even if `constr = FALSE` or an alternative constraint is specified in the 
#'   `f()` object, the same sum-to-zero constraint is applied to the model
#'   scaling.
#' 

fit_no_constr <- inla(y ~ 0 +
             f(area, model = "besag", graph = adj, hyper = hyper_area,
               scale.model = TRUE, constr = FALSE),
           data = datanull, family = "poisson",
           control.inla = list(strategy = "gaussian", int.strategy = "eb"),
           control.fixed = list(mean.intercept = 0, prec.intercept = 1),
           control.compute = list(config = TRUE))

fit_no_constr$misc$configs$constr
fit_no_constr$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]

fit_alt_constr <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                scale.model = TRUE, constr = FALSE,
                diagonal = 0,
                extraconstr = list(A = matrix(c(1, 1, 0, 0), 1), e = 3)),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

fit_alt_constr$misc$configs$constr
fit_alt_constr$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]

#' Called externally, the alternative constraint does slightly affect model 
#' scaling, and so the above confirms that the alternative constraint is not
#' used in the `scale.model` specification.

inla.scale.model(R_area, constr = list(A = matrix(c(1, 1, 1, 1), ncol = 4), e = 0))
inla.scale.model(R_area, constr = list(A = matrix(c(1, 1, 0, 0), ncol = 4), e = 3))

#' ## Default behaviour of `f(..., group = <var>)`
#' 
#' Next we review the default behaviour of the grouping option to specify product
#' smooths. The below fits a separable space-time model with an ICAR area effect
#' and RW1 time effect, otherwise using defaults.
#' 
fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                group = time, control.group = list(model = "rw1")),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))


#' ### Default constraints for grouped models
#' 
#' * By default, the a separate sum-to-zero constraint is specified for each level
#'   of the group variable.

fit$misc$configs$constr

#' * If `extraconstr=` is specified, the constraint must be the length of the
#'   number of levels for the primary variable (e.g. 4 `area`s in the example
#'   below).
#' * This constraint is repeated for each level of the `group` variable.

##+ results = "show"
area_constr <- list(A = matrix(c(1.5, 0.5, 0, 0), 1), e = 3)

fit_alt_constr <- inla(y ~ 0 +
       f(area, model = "besag", graph = adj, hyper = hyper_area,
         group = time, control.group = list(model = "rw1"),
         constr = FALSE, extraconstr = area_constr),
     data = datanull, family = "poisson",
     control.inla = list(strategy = "gaussian", int.strategy = "eb"),
     control.fixed = list(mean.intercept = 0, prec.intercept = 1),
     control.compute = list(config = TRUE))

fit_alt_constr$misc$configs$constr

#' As far as I can tell, there is no way to specify (1) a constraint for only
#' some group levels, (2) different constraints for different group levels, or
#' (3) constraints that span multiple group levels. For example, the following 
#' constraint with dimension 4 x 3 = 12 might be used to specify an overall 
#' sum-to-zero constraint (across all groups) on the space x time latent field.
#' 

wishful_constr <- list(A = matrix(1, ncol = 12), e = 0)
wishful_constr

#' But this throws an error because INLA is expecting the constraint to have 
#' dimension 4 (the number of areas) and will repeat this constraint multiple
#' times.

##+ error = TRUE
inla(y ~ 0 +
       f(area, model = "besag", graph = adj, hyper = hyper_area,
         group = time, control.group = list(model = "rw1"),
         constr = FALSE, 
         extraconstr = wishful_constr),
     data = datanull, family = "poisson")

#' More flexible constraints should be possible by specifying custom model 
#' using the `"rgeneric"` model type and manually specifying the `Cmatrix`
#' as the Kronecker product.
#' 
#' 
#' ### Model scaling with grouped models
#'  
#' * The `control.group = list(...)` default for `scale.model` is
#'   `TRUE`, even when the default for `scale.model` is `FALSE`.
#' * The main term is not scaled, consistent with `f(..., scale.model = )` default.
#'    
#'  
#'  

inla.getOption()$scale.model.default

#' Model with defaults for scale model

fit <- inla(y ~ 0 +
             f(area, model = "besag", graph = adj, hyper = hyper_area,
               group = time, control.group = list(model = "rw1")),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

R_area_scaled <- inla.scale.model(R_area, 
                                  constr = list(A = matrix(1, ncol = 4), e = 0))
R_time_scaled <- inla.scale.model(R_time, 
                                  constr = list(A = matrix(1, ncol = 3), e = 0))

round(fit$misc$configs$config[[1]]$Q[-(1:3), -(1:3)], 5)

#' This matches the Kronecker product of the unscaled area structure matrix
#' and scaled time structure matrix.

kronecker(R_time_scaled, R_area)


#' Fit both scaled `f(..., scale.model = TRUE, ..., control.group = list(..., scale.model = TRUE))`

fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                scale.model = TRUE,
                group = time,
                control.group = list(model = "rw1", scale.model = TRUE)),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

round(fit$misc$configs$config[[1]]$Q[-(1:3), -(1:3)], 5)
round(kronecker(R_time_scaled, R_area_scaled), 5)


#'
#' Fit main effect scaled and group unscaled
#' `f(..., scale.model = TRUE, ..., control.group = list(..., scale.model = FALSE))`.
#' 

fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                scale.model = TRUE,
                group = time,
                control.group = list(model = "rw1", scale.model = FALSE)),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

round(fit$misc$configs$config[[1]]$Q[-(1:3), -(1:3)], 5)
round(kronecker(R_time, R_area_scaled), 5)


#' # Confirm we can get the same thing by manually specifying `Cmatrix`.

fit <- inla(y ~ 1 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                scale.model = TRUE,
                group = time,
                control.group = list(model = "rw1", scale.model = FALSE)),
            data = data, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.compute = list(config = TRUE))


Q <- kronecker(R_time, R_area_scaled)
data$area.time <- 1:12

A <- rbind(rep(c(1, 0), c(4, 8)),
           rep(c(0, 1, 0), c(4, 4, 4)),
           rep(c(0, 1), c(8, 4)))
e <- c(0, 0, 0)


diagval <- INLA:::inla.set.f.default()$diagonal

fitQ <- inla(y ~ 1 +
               f(area.time, model = "generic0", Cmatrix = Q, hyper = hyper_area,
                 diagonal = diagval, extraconstr = list(A = A, e = e)),
             data = data, family = "poisson",
             control.inla = list(strategy = "gaussian", int.strategy = "eb"),
             control.compute = list(config = TRUE))

#' The effective parameters are equal:

fit$neffp
fitQ$neffp

#' Fixed effect and random effects are equal:

fit$summary.fixed[ , 1:2]
fitQ$summary.fixed[ , 1:2]

fit$summary.random[[1]][ , 1:3]
fitQ$summary.random[[1]][ , 1:3]

#' Marginal likelihood slightly different. Suspect difference in scaling constants?

fit$mlik
fitQ$mlik

#' # The AR1 model
#'
#' ## Basic AR1 model
#' 
#' * The raw AR1 structure matrix is scaled by $1/(1-\rho^2)$ such that the marginal
#'   variance of the precision matrix is the specified precision.

rho <- 0.8

#' The AR1 precision matrix is given by

R_ar1 <- sparseMatrix(1:5, 1:5, x = 1)
diag(R_ar1)[2:4] <- 1 + rho^2
R_ar1[cbind(1:4, 2:5)] <- -rho
R_ar1[cbind(2:5, 1:4)] <- -rho

#' For a marginal precision of 1.0, scale the raw structure matrix:

Q_ar1 <- R_ar1 * 1 / (1 - rho^2)

#' Show this matches the Q matrix constructed by INLA

hyper_ar1 <- list(prec = list(initial = log(1), fixed = TRUE),
                  rho = list(initial = log((1 + rho) / (1 - rho)), fixed = TRUE))
                  
fit <- inla(y ~ 0 +
              f(time, model = "ar1", values = 1:5, hyper = hyper_ar1),
            data = datanull[1, ], family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

fit$misc$configs$config[[1]]$Q[-1, -1]

Q_ar1

#' ## ICAR x AR1 interaction

#' Three time points 

R_ar1 <- sparseMatrix(1:3, 1:3, x = 1)
diag(R_ar1)[2] <- 1 + rho^2
R_ar1[cbind(1:2, 2:3)] <- -rho
R_ar1[cbind(2:3, 1:2)] <- -rho

Q_ar1 <- R_ar1 * 1 / (1 - rho^2)

fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                group = time,
                control.group = list(model = "ar1", hyper = hyper_ar1["rho"])),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

round(fit$misc$configs$config[[1]]$Q[-(1:3), -(1:3)], 5)
kronecker(Q_ar1, R_area)


#' ## ICAR x AR1 interaction with `scale.model = TRUE`

fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                scale.model = TRUE,
                group = time, 
                control.group = list(model = "ar1", hyper = hyper_ar1["rho"])),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

round(fit$misc$configs$config[[1]]$Q[-(1:3), -(1:3)], 5)
round(kronecker(Q_ar1, R_area_scaled), 5)



#' # `sessionInfo()`

sessionInfo()
