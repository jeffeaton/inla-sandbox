This document is to record the details of some of the internal
specifications of INLA. The purpose of this is largely to document what
INLA is doing for comparing model implementations in other software.

**Make sure to be using a version of INLA more recent than
INLA\_20.06.15 *testing* version when the internal implementation of the
RW1 and RW2 models was updated to be consistent with manual model
scaling.**

``` r
library(INLA)
```

Data
====

Create some very simple space-time data for testing INLA. Data are for
four areas in which area 2 is connected to all other areas, areas 1 and
3 are connected, and area 4 is connected only to area 2. Three time
points are simulated.

``` r
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
```

    ##   a b c d
    ## a 0 1 1 0
    ## b 1 0 1 1
    ## c 1 1 0 0
    ## d 0 1 0 0

Structure matrices for ICAR area effect and RW1 time effects.

``` r
R_area <- diag(rowSums(adj)) - adj
R_time <- t(diff(diag(3))) %*% diff(diag(3))

R_area
```

    ##    a  b  c  d
    ## a  2 -1 -1  0
    ## b -1  3 -1 -1
    ## c -1 -1  2  0
    ## d  0 -1  0  1

``` r
R_time
```

    ##      [,1] [,2] [,3]
    ## [1,]    1   -1    0
    ## [2,]   -1    2   -1
    ## [3,]    0   -1    1

Details of INLA internals
=========================

This section documents details of how certain arguments and options are
implemented by INLA such as how the small constant is added to the
diagonal, model scaling, and kronecker product for the
`f(..., group = <>)` option.

The strategy for checking is to create a ‘null’ data set with no
observations and ‘fit’ the model with fixed values for the hyper
parameters to recover the Q matrix constructed by INLA.

``` r
datanull <- data[data$area == 1, ]
datanull$y <- NA

hyper_area <- list(prec = list(initial = log(1), fixed = TRUE))
hyper_time <- list(prec = list(initial = log(1), fixed = TRUE))
```

ICAR model
----------

``` r
fit1 <- inla(y ~ 0 +
               f(area, model = "besag", graph = adj, hyper = hyper_area),
             data = datanull, family = "poisson",
             control.inla = list(strategy = "gaussian", int.strategy = "eb"),
             control.fixed = list(mean.intercept = 0, prec.intercept = 1),
             control.compute = list(config = TRUE))
```

### Extra constant added to diagonal

-   The default value for the diagonal extra constant is 1.0151131^{-5}.
    -   This is ascertained from `INLA:::inla.set.f.default()$diagonal`.
-   This is added *after* scaling the matrix by the precision parameter.

To see this, in `fit1` the fixed value for the precision is 1.0 and in
`fit2` the value for the precision is 2.0. The added value along the
diagonal in both cases is 1e-5.

``` r
hyper_area2 <- list(prec = list(initial = log(2), fixed = TRUE))

fit2 <- inla(y ~ 0 +
               f(area, model = "besag", graph = adj, hyper = hyper_area2),
             data = datanull, family = "poisson",
             control.inla = list(strategy = "gaussian", int.strategy = "eb"),
             control.fixed = list(mean.intercept = 0, prec.intercept = 1),
             control.compute = list(config = TRUE))

fit1$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]
```

    ## 4 x 4 sparse Matrix of class "dgCMatrix"
    ##                                         
    ## [1,] 488266.4 -1.00000 -1.00000  .      
    ## [2,]      .    3.00001 -1.00000 -1.00000
    ## [3,]      .    .        2.00001  .      
    ## [4,]      .    .        .        1.00001

``` r
fit2$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]
```

    ## 4 x 4 sparse Matrix of class "dgCMatrix"
    ##                                         
    ## [1,] 488268.4 -2.00000 -2.00000  .      
    ## [2,]      .    6.00001 -2.00000 -2.00000
    ## [3,]      .    .        4.00001  .      
    ## [4,]      .    .        .        2.00001

### Application of `f(..., scale.model = TRUE)`

When argument `scale.model = TRUE`, the precision matrix is scaled so
that the generalised variance is 1.

-   A sum-to-zero constraint is applied in the model scaling. (Likely a
    different constraint is applied if there are multiple connected
    components).
-   No constant is added to the diagonal before model scaling.

``` r
fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                scale.model = TRUE),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))

fit$misc$configs$constr
```

    ## $nc
    ## [1] 1
    ## 
    ## $A
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
    ## [1,]    0    0    0    1    1    1    1
    ## 
    ## $e
    ## [1] 0

``` r
fit$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]
```

    ## 4 x 4 sparse Matrix of class "dgCMatrix"
    ##                                               
    ## [1,] 488265.1 -0.3565926 -0.3565926  .        
    ## [2,]      .    1.0697879 -0.3565926 -0.3565926
    ## [3,]      .    .          0.7131953  .        
    ## [4,]      .    .          .          0.3566028

``` r
inla.scale.model(R_area, constr = list(A = matrix(1, ncol = 4), e = 0))
```

    ## 4 x 4 sparse Matrix of class "dgTMatrix"
    ##            a          b          c          d
    ## a  0.7131852 -0.3565926 -0.3565926  .        
    ## b -0.3565926  1.0697778 -0.3565926 -0.3565926
    ## c -0.3565926 -0.3565926  0.7131852  .        
    ## d  .         -0.3565926  .          0.3565926

-   Even if `constr = FALSE` or an alternative constraint is specified
    in the `f()` object, the same sum-to-zero constraint is applied to
    the model scaling.

``` r
fit_no_constr <- inla(y ~ 0 +
             f(area, model = "besag", graph = adj, hyper = hyper_area,
               scale.model = TRUE, constr = FALSE),
           data = datanull, family = "poisson",
           control.inla = list(strategy = "gaussian", int.strategy = "eb"),
           control.fixed = list(mean.intercept = 0, prec.intercept = 1),
           control.compute = list(config = TRUE))

fit_no_constr$misc$configs$constr
```

    ## NULL

``` r
fit_no_constr$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]
```

    ## 4 x 4 sparse Matrix of class "dgCMatrix"
    ##                                               
    ## [1,] 488265.1 -0.3565926 -0.3565926  .        
    ## [2,]      .    1.0697779 -0.3565926 -0.3565926
    ## [3,]      .    .          0.7131853  .        
    ## [4,]      .    .          .          0.3565926

``` r
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
```

    ## $nc
    ## [1] 1
    ## 
    ## $A
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
    ## [1,]    0    0    0    1    1    0    0
    ## 
    ## $e
    ## [1] 3

``` r
fit_alt_constr$misc$configs$config[[1]]$Q[-(1:3), -(1:3)]
```

    ## 4 x 4 sparse Matrix of class "dgCMatrix"
    ##                                               
    ## [1,] 488265.1 -0.3565926 -0.3565926  .        
    ## [2,]      .    1.0697779 -0.3565926 -0.3565926
    ## [3,]      .    .          0.7131853  .        
    ## [4,]      .    .          .          0.3565926

Called externally, the alternative constraint does slightly affect model
scaling, and so the above confirms that the alternative constraint is
not used in the `scale.model` specification.

``` r
inla.scale.model(R_area, constr = list(A = matrix(c(1, 1, 1, 1), ncol = 4), e = 0))
```

    ## 4 x 4 sparse Matrix of class "dgTMatrix"
    ##            a          b          c          d
    ## a  0.7131852 -0.3565926 -0.3565926  .        
    ## b -0.3565926  1.0697778 -0.3565926 -0.3565926
    ## c -0.3565926 -0.3565926  0.7131852  .        
    ## d  .         -0.3565926  .          0.3565926

``` r
inla.scale.model(R_area, constr = list(A = matrix(c(1, 1, 0, 0), ncol = 4), e = 3))
```

    ## 4 x 4 sparse Matrix of class "dgTMatrix"
    ##            a          b          c          d
    ## a  0.7135650 -0.3567825 -0.3567825  .        
    ## b -0.3567825  1.0703475 -0.3567825 -0.3567825
    ## c -0.3567825 -0.3567825  0.7135650  .        
    ## d  .         -0.3567825  .          0.3567825

Default behaviour of `f(..., group = <var>)`
--------------------------------------------

Next we review the default behaviour of the grouping option to specify
product smooths. The below fits a separable space-time model with an
ICAR area effect and RW1 time effect, otherwise using defaults.

``` r
fit <- inla(y ~ 0 +
              f(area, model = "besag", graph = adj, hyper = hyper_area,
                group = time, control.group = list(model = "rw1")),
            data = datanull, family = "poisson",
            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
            control.fixed = list(mean.intercept = 0, prec.intercept = 1),
            control.compute = list(config = TRUE))
```

### Default constraints for grouped models

-   By default, the a separate sum-to-zero constraint is specified for
    each level of the group variable.

``` r
fit$misc$configs$constr
```

    ## $nc
    ## [1] 3
    ## 
    ## $A
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    ## [1,]    0    0    0    1    1    1    1    0    0     0     0     0     0     0
    ## [2,]    0    0    0    0    0    0    0    1    1     1     1     0     0     0
    ## [3,]    0    0    0    0    0    0    0    0    0     0     0     1     1     1
    ##      [,15]
    ## [1,]     0
    ## [2,]     0
    ## [3,]     1
    ## 
    ## $e
    ## [1] 0 0 0

-   If `extraconstr=` is specified, the constraint must be the length of
    the number of levels for the primary variable (e.g. 4 `area`s in the
    example below).
-   This constraint is repeated for each level of the `group` variable.

``` r
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
```

    ## $nc
    ## [1] 3
    ## 
    ## $A
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    ## [1,]    0    0    0  1.5  0.5    0    0  0.0  0.0     0     0   0.0   0.0     0
    ## [2,]    0    0    0  0.0  0.0    0    0  1.5  0.5     0     0   0.0   0.0     0
    ## [3,]    0    0    0  0.0  0.0    0    0  0.0  0.0     0     0   1.5   0.5     0
    ##      [,15]
    ## [1,]     0
    ## [2,]     0
    ## [3,]     0
    ## 
    ## $e
    ## [1] 3 3 3

As far as I can tell, there is no way to specify (1) a constraint for
only some group levels, (2) different constraints for different group
levels, or (3) constraints that span multiple group levels. For example,
the following constraint with dimension 4 x 3 = 12 might be used to
specify an overall sum-to-zero constraint (across all groups) on the
space x time latent field.

``` r
wishful_constr <- list(A = matrix(1, ncol = 12), e = 0)
wishful_constr
```

    ## $A
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    1    1    1    1    1    1    1    1    1     1     1     1
    ## 
    ## $e
    ## [1] 0

But this throws an error because INLA is expecting the constraint to
have dimension 4 (the number of areas) and will repeat this constraint
multiple times.

``` r
inla(y ~ 0 +
       f(area, model = "besag", graph = adj, hyper = hyper_area,
         group = time, control.group = list(model = "rw1"),
         constr = FALSE, 
         extraconstr = wishful_constr),
     data = datanull, family = "poisson")
```

    ## Error in inla(y ~ 0 + f(area, model = "besag", graph = adj, hyper = hyper_area, : 
    ##  ncol in matrix A(extraconstr) does not correspont to the length of f: 12 4

More flexible constraints should be possible by specifying custom model
using the `"rgeneric"` model type and manually specifying the `Cmatrix`
as the Kronecker product.

### Model scaling with grouped models

-   The `control.group = list(...)` default for `scale.model` is `TRUE`,
    even when the default for `scale.model` is `FALSE`.
-   The main term is not scaled, consistent with
    `f(..., scale.model = )` default.

``` r
inla.getOption()$scale.model.default
```

    ## [1] FALSE

Model with defaults for scale model

``` r
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
```

    ## 12 x 12 sparse Matrix of class "dgCMatrix"
    ##                                                                         
    ##  [1,] 162755.6 -0.40934 -0.40934  .           -0.81867  0.40934  0.40934
    ##  [2,]      .    1.22802 -0.40934 -0.40934      0.40934 -1.22801  0.40934
    ##  [3,]      .    .        0.81868  .            0.40934  0.40934 -0.81867
    ##  [4,]      .    .        .        0.40935      .        0.40934  .      
    ##  [5,]      .    .        .        .       162756.42878 -0.81867 -0.81867
    ##  [6,]      .    .        .        .            .        2.45603 -0.81867
    ##  [7,]      .    .        .        .            .        .        1.63736
    ##  [8,]      .    .        .        .            .        .        .      
    ##  [9,]      .    .        .        .            .        .        .      
    ## [10,]      .    .        .        .            .        .        .      
    ## [11,]      .    .        .        .            .        .        .      
    ## [12,]      .    .        .        .            .        .        .      
    ##                                                       
    ##  [1,]  .            .        .        .        .      
    ##  [2,]  0.40934      .        .        .        .      
    ##  [3,]  .            .        .        .        .      
    ##  [4,] -0.40934      .        .        .        .      
    ##  [5,]  .           -0.81867  0.40934  0.40934  .      
    ##  [6,] -0.81867      0.40934 -1.22801  0.40934  0.40934
    ##  [7,]  .            0.40934  0.40934 -0.81867  .      
    ##  [8,]  0.81868      .        0.40934  .       -0.40934
    ##  [9,]  .       162755.61010 -0.40934 -0.40934  .      
    ## [10,]  .            .        1.22802 -0.40934 -0.40934
    ## [11,]  .            .        .        0.81868  .      
    ## [12,]  .            .        .        .        0.40935

This matches the Kronecker product of the unscaled area structure matrix
and scaled time structure matrix.

``` r
kronecker(R_time_scaled, R_area)
```

    ## 12 x 12 sparse Matrix of class "dgTMatrix"
    ##                                                                        
    ##  [1,]  0.8186736 -0.4093368 -0.4093368  .         -0.8186736  0.4093368
    ##  [2,] -0.4093368  1.2280105 -0.4093368 -0.4093368  0.4093368 -1.2280105
    ##  [3,] -0.4093368 -0.4093368  0.8186736  .          0.4093368  0.4093368
    ##  [4,]  .         -0.4093368  .          0.4093368  .          0.4093368
    ##  [5,] -0.8186736  0.4093368  0.4093368  .          1.6373473 -0.8186736
    ##  [6,]  0.4093368 -1.2280105  0.4093368  0.4093368 -0.8186736  2.4560209
    ##  [7,]  0.4093368  0.4093368 -0.8186736  .         -0.8186736 -0.8186736
    ##  [8,]  .          0.4093368  .         -0.4093368  .         -0.8186736
    ##  [9,]  .          .          .          .         -0.8186736  0.4093368
    ## [10,]  .          .          .          .          0.4093368 -1.2280105
    ## [11,]  .          .          .          .          0.4093368  0.4093368
    ## [12,]  .          .          .          .          .          0.4093368
    ##                                                                        
    ##  [1,]  0.4093368  .          .          .          .          .        
    ##  [2,]  0.4093368  0.4093368  .          .          .          .        
    ##  [3,] -0.8186736  .          .          .          .          .        
    ##  [4,]  .         -0.4093368  .          .          .          .        
    ##  [5,] -0.8186736  .         -0.8186736  0.4093368  0.4093368  .        
    ##  [6,] -0.8186736 -0.8186736  0.4093368 -1.2280105  0.4093368  0.4093368
    ##  [7,]  1.6373473  .          0.4093368  0.4093368 -0.8186736  .        
    ##  [8,]  .          0.8186736  .          0.4093368  .         -0.4093368
    ##  [9,]  0.4093368  .          0.8186736 -0.4093368 -0.4093368  .        
    ## [10,]  0.4093368  0.4093368 -0.4093368  1.2280105 -0.4093368 -0.4093368
    ## [11,] -0.8186736  .         -0.4093368 -0.4093368  0.8186736  .        
    ## [12,]  .         -0.4093368  .         -0.4093368  .          0.4093368

Fit both scaled
`f(..., scale.model = TRUE, ..., control.group = list(..., scale.model = TRUE))`

``` r
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
```

    ## 12 x 12 sparse Matrix of class "dgCMatrix"
    ##                                                                         
    ##  [1,] 162755.1 -0.14597 -0.14597  .           -0.29193  0.14597  0.14597
    ##  [2,]      .    0.43791 -0.14597 -0.14597      0.14597 -0.43790  0.14597
    ##  [3,]      .    .        0.29194  .            0.14597  0.14597 -0.29193
    ##  [4,]      .    .        .        0.14598      .        0.14597  .      
    ##  [5,]      .    .        .        .       162755.37530 -0.29193 -0.29193
    ##  [6,]      .    .        .        .            .        0.87581 -0.29193
    ##  [7,]      .    .        .        .            .        .        0.58388
    ##  [8,]      .    .        .        .            .        .        .      
    ##  [9,]      .    .        .        .            .        .        .      
    ## [10,]      .    .        .        .            .        .        .      
    ## [11,]      .    .        .        .            .        .        .      
    ## [12,]      .    .        .        .            .        .        .      
    ##                                                       
    ##  [1,]  .            .        .        .        .      
    ##  [2,]  0.14597      .        .        .        .      
    ##  [3,]  .            .        .        .        .      
    ##  [4,] -0.14597      .        .        .        .      
    ##  [5,]  .           -0.29193  0.14597  0.14597  .      
    ##  [6,] -0.29193      0.14597 -0.43790  0.14597  0.14597
    ##  [7,]  .            0.14597  0.14597 -0.29193  .      
    ##  [8,]  0.29194      .        0.14597  .       -0.14597
    ##  [9,]  .       162755.08336 -0.14597 -0.14597  .      
    ## [10,]  .            .        0.43791 -0.14597 -0.14597
    ## [11,]  .            .        .        0.29194  .      
    ## [12,]  .            .        .        .        0.14598

``` r
round(kronecker(R_time_scaled, R_area_scaled), 5)
```

    ## 12 x 12 sparse Matrix of class "dgTMatrix"
    ##                                                                              
    ##  [1,]  0.29193 -0.14597 -0.14597  .       -0.29193  0.14597  0.14597  .      
    ##  [2,] -0.14597  0.43790 -0.14597 -0.14597  0.14597 -0.43790  0.14597  0.14597
    ##  [3,] -0.14597 -0.14597  0.29193  .        0.14597  0.14597 -0.29193  .      
    ##  [4,]  .       -0.14597  .        0.14597  .        0.14597  .       -0.14597
    ##  [5,] -0.29193  0.14597  0.14597  .        0.58387 -0.29193 -0.29193  .      
    ##  [6,]  0.14597 -0.43790  0.14597  0.14597 -0.29193  0.87580 -0.29193 -0.29193
    ##  [7,]  0.14597  0.14597 -0.29193  .       -0.29193 -0.29193  0.58387  .      
    ##  [8,]  .        0.14597  .       -0.14597  .       -0.29193  .        0.29193
    ##  [9,]  .        .        .        .       -0.29193  0.14597  0.14597  .      
    ## [10,]  .        .        .        .        0.14597 -0.43790  0.14597  0.14597
    ## [11,]  .        .        .        .        0.14597  0.14597 -0.29193  .      
    ## [12,]  .        .        .        .        .        0.14597  .       -0.14597
    ##                                          
    ##  [1,]  .        .        .        .      
    ##  [2,]  .        .        .        .      
    ##  [3,]  .        .        .        .      
    ##  [4,]  .        .        .        .      
    ##  [5,] -0.29193  0.14597  0.14597  .      
    ##  [6,]  0.14597 -0.43790  0.14597  0.14597
    ##  [7,]  0.14597  0.14597 -0.29193  .      
    ##  [8,]  .        0.14597  .       -0.14597
    ##  [9,]  0.29193 -0.14597 -0.14597  .      
    ## [10,] -0.14597  0.43790 -0.14597 -0.14597
    ## [11,] -0.14597 -0.14597  0.29193  .      
    ## [12,]  .       -0.14597  .        0.14597

Fit main effect scaled and group unscaled
`f(..., scale.model = TRUE, ..., control.group = list(..., scale.model = FALSE))`.

``` r
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
```

    ## 12 x 12 sparse Matrix of class "dgCMatrix"
    ##                                                                         
    ##  [1,] 162755.5 -0.35659 -0.35659  .           -0.71319  0.35659  0.35659
    ##  [2,]      .    1.06979 -0.35659 -0.35659      0.35659 -1.06978  0.35659
    ##  [3,]      .    .        0.71320  .            0.35659  0.35659 -0.71319
    ##  [4,]      .    .        .        0.35660      .        0.35659  .      
    ##  [5,]      .    .        .        .       162756.21780 -0.71319 -0.71319
    ##  [6,]      .    .        .        .            .        2.13957 -0.71319
    ##  [7,]      .    .        .        .            .        .        1.42638
    ##  [8,]      .    .        .        .            .        .        .      
    ##  [9,]      .    .        .        .            .        .        .      
    ## [10,]      .    .        .        .            .        .        .      
    ## [11,]      .    .        .        .            .        .        .      
    ## [12,]      .    .        .        .            .        .        .      
    ##                                                       
    ##  [1,]  .            .        .        .        .      
    ##  [2,]  0.35659      .        .        .        .      
    ##  [3,]  .            .        .        .        .      
    ##  [4,] -0.35659      .        .        .        .      
    ##  [5,]  .           -0.71319  0.35659  0.35659  .      
    ##  [6,] -0.71319      0.35659 -1.06978  0.35659  0.35659
    ##  [7,]  .            0.35659  0.35659 -0.71319  .      
    ##  [8,]  0.71320      .        0.35659  .       -0.35659
    ##  [9,]  .       162755.50461 -0.35659 -0.35659  .      
    ## [10,]  .            .        1.06979 -0.35659 -0.35659
    ## [11,]  .            .        .        0.71320  .      
    ## [12,]  .            .        .        .        0.35660

``` r
round(kronecker(R_time, R_area_scaled), 5)
```

    ## 12 x 12 sparse Matrix of class "dgTMatrix"
    ##                                                                              
    ##  [1,]  0.71319 -0.35659 -0.35659  .       -0.71319  0.35659  0.35659  .      
    ##  [2,] -0.35659  1.06978 -0.35659 -0.35659  0.35659 -1.06978  0.35659  0.35659
    ##  [3,] -0.35659 -0.35659  0.71319  .        0.35659  0.35659 -0.71319  .      
    ##  [4,]  .       -0.35659  .        0.35659  .        0.35659  .       -0.35659
    ##  [5,] -0.71319  0.35659  0.35659  .        1.42637 -0.71319 -0.71319  .      
    ##  [6,]  0.35659 -1.06978  0.35659  0.35659 -0.71319  2.13956 -0.71319 -0.71319
    ##  [7,]  0.35659  0.35659 -0.71319  .       -0.71319 -0.71319  1.42637  .      
    ##  [8,]  .        0.35659  .       -0.35659  .       -0.71319  .        0.71319
    ##  [9,]  .        .        .        .       -0.71319  0.35659  0.35659  .      
    ## [10,]  .        .        .        .        0.35659 -1.06978  0.35659  0.35659
    ## [11,]  .        .        .        .        0.35659  0.35659 -0.71319  .      
    ## [12,]  .        .        .        .        .        0.35659  .       -0.35659
    ##                                          
    ##  [1,]  .        .        .        .      
    ##  [2,]  .        .        .        .      
    ##  [3,]  .        .        .        .      
    ##  [4,]  .        .        .        .      
    ##  [5,] -0.71319  0.35659  0.35659  .      
    ##  [6,]  0.35659 -1.06978  0.35659  0.35659
    ##  [7,]  0.35659  0.35659 -0.71319  .      
    ##  [8,]  .        0.35659  .       -0.35659
    ##  [9,]  0.71319 -0.35659 -0.35659  .      
    ## [10,] -0.35659  1.06978 -0.35659 -0.35659
    ## [11,] -0.35659 -0.35659  0.71319  .      
    ## [12,]  .       -0.35659  .        0.35659

Confirm we can get the same thing by manually specifying `Cmatrix`.
===================================================================

``` r
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
```

The effective parameters are equal:

``` r
fit$neffp
```

    ##                                       [,1]
    ## Expectected  number of parameters 7.576113
    ## Stdev of the number of parameters 0.000000
    ## Number of equivalent replicates   1.583926

``` r
fitQ$neffp
```

    ##                                       [,1]
    ## Expectected  number of parameters 7.576115
    ## Stdev of the number of parameters 0.000000
    ## Number of equivalent replicates   1.583925

Fixed effect and random effects are equal:

``` r
fit$summary.fixed[ , 1:2]
```

    ##                  mean        sd
    ## (Intercept) 0.8251967 0.1972447

``` r
fitQ$summary.fixed[ , 1:2]
```

    ##                  mean        sd
    ## (Intercept) 0.8251967 0.1972447

``` r
fit$summary.random[[1]][ , 1:3]
```

    ##    ID        mean        sd
    ## 1   1 -0.94219909 0.6243636
    ## 2   2 -0.00137828 0.4842483
    ## 3   3  0.24255352 0.4697330
    ## 4   4  0.70102386 0.4345768
    ## 5   1 -0.75627048 0.5427202
    ## 6   2  0.25072429 0.4145179
    ## 7   3  0.44648383 0.4143926
    ## 8   4  0.05906236 0.4967552
    ## 9   1  0.28635644 0.4684865
    ## 10  2 -0.24403882 0.5025126
    ## 11  3  0.10077233 0.4900527
    ## 12  4 -0.14308995 0.5497946

``` r
fitQ$summary.random[[1]][ , 1:3]
```

    ##    ID        mean        sd
    ## 1   1 -0.94219909 0.6243636
    ## 2   2 -0.00137828 0.4842483
    ## 3   3  0.24255352 0.4697330
    ## 4   4  0.70102386 0.4345768
    ## 5   5 -0.75627048 0.5427202
    ## 6   6  0.25072429 0.4145179
    ## 7   7  0.44648383 0.4143926
    ## 8   8  0.05906236 0.4967552
    ## 9   9  0.28635644 0.4684865
    ## 10 10 -0.24403882 0.5025126
    ## 11 11  0.10077233 0.4900527
    ## 12 12 -0.14308995 0.5497947

Marginal likelihood slightly different. Suspect difference in scaling
constants?

``` r
fit$mlik
```

    ##                                            [,1]
    ## log marginal-likelihood (integration) -23.10433
    ## log marginal-likelihood (Gaussian)    -23.10433

``` r
fitQ$mlik
```

    ##                                            [,1]
    ## log marginal-likelihood (integration) -25.86114
    ## log marginal-likelihood (Gaussian)    -25.86114

`sessionInfo()`
===============

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Catalina 10.15.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] INLA_20.06.15 foreach_1.5.0 sp_1.4-2      Matrix_1.2-18
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] codetools_0.2-16   lattice_0.20-41    digest_0.6.25      grid_3.6.3        
    ##  [5] MatrixModels_0.4-1 magrittr_1.5       evaluate_0.14      highr_0.8         
    ##  [9] rlang_0.4.6        stringi_1.4.6      rmarkdown_2.2      splines_3.6.3     
    ## [13] iterators_1.0.12   tools_3.6.3        stringr_1.4.0      xfun_0.14         
    ## [17] yaml_2.2.1         compiler_3.6.3     htmltools_0.5.0    knitr_1.28
