rm(list = ls())
require(EQL)
require(data.table)
require(ggplot2)

setwd('~/git/abbg/R')
source('fun.model_solver_nl_kg.r')
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

#grid dimension
p$nbin  = 50    #permanent component
p$neps  = 49     #transitory component
p$ngpa  = 50     #asset
p$ngpm  = 49     #average earnings points
p$ngpp  = p$ngpm * p$nbin * p$neps      #pension points

#asset     
p$amax = 300000   
p$pexpgrid = 0.18  #approaches linear as goes to 0, approaches L shaped as goes to Inf
              
#utility
p$R       = 1.06             #(* Gross interest rate *)
p$beta    = 0.93             #(* Discount factor *)
p$rho     = 2                #(* Coefficient of Relative Risk Aversion *)
p$cfloor = -100000000000

#sim
p$nsim  = 50000     # Number of people to simulate
p$N     = 999999

#age
p$age_min = 30
p$age_re  = 66
p$age_max = 90
p$nage  = (p$age_re - p$age_min)/2
p$Tret   = (p$age_max - p$age_re)/2 +1 #periods after retirement

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