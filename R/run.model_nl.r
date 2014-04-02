rm(list = ls())
require(EQL)
require(data.table)
require(ggplot2)

setwd('~/git/abbg/R')
source('fun.model_solver_sim.r')
set.seed(123)

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
p$meanAGE    = meanAGE    
p$stdAGE     = stdAGE
p$Vectau     = Vectau
p$Ntau       = Ntau
p$meanY      = meanY
p$stdY       = stdY

detach(data)

p$nage  = (60-30)/2 + 1
p$nbin  = 23

#compute the transition matrix of persistent income
eta <- comp.eta(p)


p$rho     = 2                #(* Coefficient of Relative Risk Aversion *)
p$R       = 1.03             #(* Gross interest rate *)
p$beta    = 0.96             #(* Discount factor *)
p$nP      = 6                #(* Permanent shock Number of points in the discrete approximation to lognormal dist *)
p$nT      = 6                #(* Transitory shock Number of points in the discrete approximation to lognormal dist *)
p$sigP    = 0.1              #(* Permanent shock Standard deviation of lognormal distribution *)
p$sigT    = 0.2              #(* Transitory shock Standard deviation of lognormal distribution *)
#p$pUnemp  = 0 #0.5/100      #Probability of unemployment (when unemployed inc level is zero) 
p$aMin  = 1e-6            #(* Lower bound for GothicAVec *)
p$aMax  = 30            #(* Maximum point in GothicAVec *)
p$aHuge = 9000 
p$n = 50                    #(* Number of points in GothicAVec *)
p$NumOfPeople  = 10000     # Number of people to simulate
p$NumOfPeriodsToSimulate = 41   #Length of life in simulation (simulate until age 60)

start_time = proc.time()[3]  
model  <- comp.solveModel(p)
cat(paste('\ntotal seconds to compute Cons rule: ' , proc.time()[3] -  start_time ))

start_time = proc.time()[3]  
moments <- comp.moments(p, model) 
cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))


