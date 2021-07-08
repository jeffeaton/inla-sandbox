library(tidyverse)
library(INLA)
library(TMB)
# precompile()

data(Munich)
g = system.file("demodata/munich.graph", package="INLA")

# INLA
# -----------------------------------------------------------------------------

Q = inla.graph2matrix(g)

formula = rent ~ f(location,
  model="besag",
  graph.file=Q,
  constr=TRUE,
  scale.model=FALSE, # still inla.scale.model!
  adjust.for.con.comp=FALSE, # sum over whole graph
  hyper=list(prec=list(prior="pc.prec", param=c(1, 0.01))) # match tmb
)
mod = inla(formula,data=Munich)

# TMB.
# -----------------------------------------------------------------------------

# dyn.unload(dynlib('icar'))
compile('icar.cpp')
dyn.load(dynlib('icar'))
invisible(config(tape.parallel=FALSE, DLL='icar'))

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
    log_icar_tau=log(0.5^-2),
    log_sd_gauss=log(2)
  ),
  random='icar_vec', 
  silent=0, 
  DLL = 'icar'
)

fit <- nlminb(munich$par, munich$fn, munich$gr)
rp <- sdreport(munich, fit$par)
re <- summary.sdreport(rp, 'random')

summary(rp, 'fixed')[2, ]
mod$summary.fixed[, 1:2]

summary.sdreport(rp, 'report')
mod$summary.hyperpar[,1:2]

plotl(re[,1])
lines(mod$summary.random$location$mean, col=2)

plot(re[,1], mod$summary.random$location$mean, asp=1, xlab='tmb', ylab='inla')
abline(0, 1)
