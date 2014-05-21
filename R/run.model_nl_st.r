rm(list = ls())
require(EQL)
require(data.table)
require(ggplot2)

setwd('~/git/abbg/R')
source('fun.model_solver_nl_st.r')

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
p$nbin  = 100
p$neps  = 23

#sim
p$nsim  = 999     # Number of people to simulate
p$N     = 999999

#age
#p$ntr   = 8  #periods after retirement
for (ai in 30:55){
  p$age_min = ai
  p$age_max = ai+10
  p$nage  = 6

  set.seed(77)
  #start_time = proc.time()[3]  
  #model  <- comp.solveModel(p)
  #cat(paste('\ntotal seconds to compute Cons rule: ' , proc.time()[3] -  start_time ))

  #start_time = proc.time()[3]  
  #moments <- comp.moments(p, model) 
  #cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))

  start_time = proc.time()[3]  
  etaeps  <- comp.income(p)
  cat(paste('\ntotal seconds to compute income: ' , proc.time()[3] -  start_time ))

  savename <- paste('cohort',p$age_min,'.dat',sep='')
  save(etaeps,p,file=savename)  

}
#source('fun.sim.data.r')
#sim <- sim.small.sample(p$age_min)

#sim <- sim.origin.sample()
#sim[,lc:=log(consumption)]  #log consumption
#sim[,cage:=as.factor(age) ] #create age dummy
#sim[,lcr:=lm(lc~cage)$residuals]
#sim$cage <- NULL
#simdata <-data.matrix(sim)
#save(simdata,file='simdata.dat')
#require(R.matlab)
#writeMat('simdata.mat',simdata=simdata)
#require(plyr)
#medeta <- ddply(sim, ~age,summarise,medeta=median(eta))

#persis <- sim.persis(p,sim)
