# Simulate the life cycle model with nonlinear rule.

rm(list = ls())
require(EQL)
require(data.table)
require(ggplot2)

#setwd('~/git/abbg/R')
source('~/git/abbg/R3/fun.model.solver.nl.r')
set.seed(77)

# SETTIG PARAMETERS
p <- list()

load('~/git/abbg/R/mat.dat')
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
p$ngpe  = 99#49     #transitory component
p$ngpz  = 100#50    #permanent component
p$ngpa  = 50     #asset
p$ngpm  = 100#50     #average earnings points

#age
p$age_min = 30
p$age_re  = 66
p$age_max = 90
p$Twork  = (p$age_re - p$age_min)/2
p$Tret   = (p$age_max - p$age_re)/2 +1 #periods after retirement

#asset
p$amax  = 3e5
p$aHuge = 3e7
p$pexpgrid = 0.18  #approaches linear as goes to 0, approaches L shaped as goes to Inf

#sim
p$nsim  = 50000     # Number of people to simulate
p$N     = 999999

#utility
p$R       = 1.06             #(* Gross interest rate *)
p$bet    = 0.93             #(* Discount factor *)
p$gam     = 2                #(* Coefficient of Relative Risk Aversion *)

#BORROWING LIMIT: SET TO VERY LARGE NEGATIVE NO. FOR NBL
p$borrowlim = 0 #0.0 #-100000000.0

start_time = proc.time()[3]
model  <- comp.solveModel(p)
cat(paste('\ntotal seconds to compute Cons rule: ' , proc.time()[3] -  start_time ))

save(model,file='only_model_100gp.dat')

#load('only_model2.dat')
start_time = proc.time()[3]
moments <- comp.moments(p, model)
cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))
