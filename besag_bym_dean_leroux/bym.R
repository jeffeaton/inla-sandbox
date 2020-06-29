# The purpose is checking whether the implementation of BYM model in TMB
# yields the same estimate of latent effects

library(tidyverse)
library(INLA)
library(TMB)
# precompile()

data(Munich)
g = system.file("demodata/munich.graph", package="INLA")

# INLA
# -----------------------------------------------------------------------------

Q = inla.graph2matrix(g)
n = nrow(Q)

formula = rent ~ f(location,
  model="bym",
  graph.file=Q,
  constr=TRUE,
  scale.model=TRUE, # still inla.scale.model!
  adjust.for.con.comp=FALSE, # sum over whole graph
  hyper=list(
    theta1=list(prior="pc.prec", param=c(1, 0.01)),
    theta2=list(prior="pc.prec", param=c(1, 0.01)) # huge precision with default
                                                   # prior, practically ignoring
                                                   # spatial effect
  ) # match tmb
)
mod = inla(formula,data=Munich)

# TMB.
# -----------------------------------------------------------------------------

# dyn.unload(dynlib('icar'))
compile('bym.cpp')
dyn.load(dynlib('bym'))
invisible(config(tape.parallel=FALSE, DLL='bym'))

Q = inla.graph2matrix(g)
nnbs = rowSums(Q)-1
Q[Q==1] = -1
diag(Q) = nnbs

# INLA scales this even if scale.model=FALSE?
Q = inla.scale.model(Q, constr = list(A = matrix(1, 1, nrow(Q)), e=0))

munich = MakeADFun(
  data=list(
    rent=Munich$rent,
    location=Munich$location-1,
    R_=as.matrix(Q)
  ),
  parameter=list(
    intercept=1,
    icar_vec=rnorm(nrow(Q)),
    iid_vec=rnorm(nrow(Q)),
    log_icar_tau=log(0.5^-2),
    log_iid_tau=log(0.5^-2),
    log_sd_gauss=log(2)
  ),
  random=c('icar_vec', 'iid_vec'), 
  silent=0, 
  DLL = 'bym'
)

fit <- nlminb(munich$par, munich$fn, munich$gr)
fit
rp <- sdreport(munich, fit$par)
re <- summary.sdreport(rp, 'random')

# Fixed effects OK
rbind(
  tmb=summary(rp, 'fixed')['intercept', ], 
  inla=mod$summary.fixed[, 1:2]
)

# Random effects
# Different precision but...
summary.sdreport(rp, 'report')
mod$summary.hyperpar[,1:2]
# ...the same sum, should be OK as INLA BYM estimate u+v and u get joint u/v
# marginals, would be nice to replicate
uv_inla <- mod$summary.random$location$mean[1:nrow(Q)]
uv_tmb  <- re[1:nrow(Q),1]+re[nrow(Q)+1:nrow(Q),1]

# sum effects
par(mfrow=c(1,2))
plot(uv_tmb, col='red', type='l')
lines(uv_inla)
plot(uv_tmb, uv_inla)
abline(0, 1)

# spatial only
u_inla <- mod$summary.random$location$mean[nrow(Q)+1:nrow(Q)]
u_tmb  <- re[nrow(Q)+1:nrow(Q),1]
plot(uv_tmb, col='red', type='l')
lines(uv_inla)
plot(uv_tmb, uv_inla)
abline(0, 1)

system('rm *.*o')