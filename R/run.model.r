rm(list = ls())
setwd('~/git/abbg/R')
source('fun.model_solver_sim.r')

# SETTIG PARAMETERS
p <- list()

p$PeriodsToSolve  = 90-25

p$rho     = 2                #(* Coefficient of Relative Risk Aversion *)
p$Rhat       = 1.03             #(* Gross interest rate *)
p$beta    = 1/p$R            #(* Discount factor *)
p$nP      = 6                #(* Permanent shock Number of points in the discrete approximation to lognormal dist *)
p$nT      = 6                #(* Transitory shock Number of points in the discrete approximation to lognormal dist *)
p$sigP    = 0.128              #(* Permanent shock Standard deviation of lognormal distribution *)
p$sigT    = 0.164              #(* Transitory shock Standard deviation of lognormal distribution *)

p$aMin  = 1e-5            #(* Lower bound for GothicAVec *)
p$aMax  = 4            #(* Maximum point in GothicAVec *)
p$aHuge = 9000 
p$n = 20                    #(* Number of points in GothicAVec *)
