rm(list = ls())
setwd('~/git/abbg/R')
source('fun.model_solver_sim.r')

# SETTIG PARAMETERS
p <- list()

p$PeriodsToSolve  = 90-25

p$rho     = 2                #(* Coefficient of Relative Risk Aversion *)
p$Rhat    = 1.03             #(* Gross interest rate *)
p$beta    = 0.96             #(* Discount factor *)
p$nP      = 6                #(* Permanent shock Number of points in the discrete approximation to lognormal dist *)
p$nT      = 6                #(* Transitory shock Number of points in the discrete approximation to lognormal dist *)
p$sigP    = 0.1              #(* Permanent shock Standard deviation of lognormal distribution *)
p$sigT    = 0.2              #(* Transitory shock Standard deviation of lognormal distribution *)
#p$pUnemp  = 0 #0.5/100      #Probability of unemployment (when unemployed inc level is zero) 
p$aMin  = 1e-6            #(* Lower bound for GothicAVec *)
p$aMax  = 4            #(* Maximum point in GothicAVec *)
p$aHuge = 9000 
p$n = 20                    #(* Number of points in GothicAVec *)
p$NumOfPeople  = 1000     # Number of people to simulate
p$NumOfPeriodsToSimulate = 40   #Length of life in simulation (simulate until age 60)

start_time = proc.time()[3]  
model  <- comp.solveModel(p)
cat(paste('\ntotal seconds to compute Cons rule: ' , proc.time()[3] -  start_time ))

start_time = proc.time()[3]  
moments <- comp.moments(p, model) 
cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))


