# 
# 
# This is an implementation of Dean modification of BYM 
# Dean, C. B., M. D. Ugarte, and A. F. Militino. ‘Detecting Interaction Between Random Region and Fixed Age Effects in Disease Mapping’. Biometrics 57, no. 1 (2001): 197–202. https://doi.org/10.1111/j.0006-341X.2001.00197.x.
# 
#  

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

# dyn.unload(dynlib('bym'))
compile('dean.cpp')
dyn.load(dynlib('dean'))
invisible(config(tape.parallel=FALSE, DLL='dean'))

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
    log_tau=log(0.5^-2),
    logit_phi=logit(.5),
    log_sd_gauss=log(2)
  ),
  random=c('icar_vec', 'iid_vec'), 
  silent=0, 
  DLL = 'dean'
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
srp = summary.sdreport(rp, 'report')
srp[c('t_gau', 't_icar', 't_iid', 'tau', 'phi'), ]
mod$summary.hyperpar[,1:2]

uv_inla <- mod$summary.random$location$mean[1:nrow(Q)]
uv_tmb  <- srp[grep('^uv$', rownames(srp)), 1]

# sum effects
par(mfrow=c(1,2))
plot(uv_tmb, col='red', type='l')
lines(uv_inla)
plot(uv_tmb, uv_inla)
abline(0, 1)

# spatial only
u_inla <- mod$summary.random$location$mean[nrow(Q)+1:nrow(Q)]
u_tmb  <- re[1:nrow(Q),1]
plot(u_tmb, col='red', type='l')
lines(u_inla)
plot(u_tmb, u_inla)
abline(0, 1)

