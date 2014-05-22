rm(list = ls())
require(EQL)
require(data.table)
require(ggplot2)

#setwd('~/git/abbg/R')
source('fun.model_solver_nl_st.r')
set.seed(77)

# SETTIG PARAMETERS
p <- list()

load('mat.dat')
attach(data)

p$Resqtrue = Resqfinal
p$Resqtrue_e0 = Resqfinal.e0
p$Resqtrue_eps = Resqfinal.eps
p$b1true = b1
p$bLtrue = bL
p$b1true_e0 = b1.e0
p$bLtrue_e0 = bL.e0
p$b1true_eps = b1.eps       
p$bLtrue_eps = bL.eps

p$K1         = K1
p$K2         = K2
p$K3         = K3
p$K4         = K4
p$meanAGE    = meanAGE    
p$stdAGE     = stdAGE
p$Vectau     = Vectau
p$Ntau       = Ntau
p$meanY      = meanY
p$stdY       = stdY
p$T          = T

detach(data)

#income node
p$nbin  = 200
p$neps  = 199

#sim
p$nsim  = 999     # Number of people to simulate
p$N     = 999999

p$firstiniage <- 30
p$lastiniage <- 55
p$nage  = 6


#sim <- sim.origin.sample()
#sim[,lc:=log(consumption)]  #log consumption
#sim[,cage:=as.factor(age) ] #create age dummy
#sim[,lcr:=lm(lc~cage)$residuals]
#sim$cage <- NULL
#require(plyr)
#medeta <- ddply(sim, ~age,summarise,medeta=median(eta))

#simdata <-data.matrix(sim)
#save(simdata,file='simdata.dat')
#require(R.matlab)
#writeMat('simdata.mat',simdata=simdata)

#sim$t <-1:p$T
#sim_y <- sim[,c("id","t","Y"),with=F]
#wide_y <- reshape(sim_y, idvar='id', timevar='t', direction='wide') #each person has all age in a row
#Y <-data.matrix(wide_y[,-1,with=F])
#writeMat('Y.mat',Y=Y)


#persis <- sim.persis(p,sim)
