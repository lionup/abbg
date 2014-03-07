rm(list = ls())
setwd('~/git/abbg/R')
source('fun.model_solver_sim.r')

# SETTIG PARAMETERS
p <- list()

p$PeriodsToSolve = 21

p$Rho     = 2                #(* Coefficient of Relative Risk Aversion *)
p$Beta    = 0.96             #(* Discount factor *)
p$NumOfThetaShockPoints = 7  #(* Number of points in the discrete approximation to lognormal dist *)
p$Sigma   = 0.2              #(* Standard deviation of lognormal distribution *)
p$RFree   = 1.02             #(* Gross interest rate *)

p$Gamma   = 1.00             # (* Permanent income growth factor *)

p$mMin    = 0                #(* Minimum point in mVec *)
p$mMax    = 4                #(* Maximum point in mVec *)
p$mHuge   = 5                #mHuge is a point so that extrapolation is not needed
p$NumOfmPts = 5              #(* Number of points in mVec *) %Perhaps change back to 5

p$GothicAMin  = 0            #(* Lower bound for GothicAVec *)
p$GothicAMax  = 4            #(* Maximum point in GothicAVec *)
p$GothicAHuge = 9000 
p$NumOfGothicAPts = 5        #(* Number of points in GothicAVec *)

p$PeriodsSolved = 0          #(* Number of periods back from T for which the model has been solved *)

p$Constrained = 0            # Constrained

p$MC = 0                     # Indicates if if is the multicontrol problem
