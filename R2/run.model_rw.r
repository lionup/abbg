#cononical random walk model
#Initial permanent component: sig2eta1 =  0.15
#Permanent shock sig2v =   0.01
#Transitory shock: sig2eps =   0.05

rm(list = ls())
setwd('~/git/abbg/R2')
source('fun.model_solver_rw.r')

# SETTIG PARAMETERS
p <- list()

#
p$age_min = 25
p$age_re  = p$age_min+35 #first period income drop at age 60 
p$age_max = p$age_min+69 #last period 94
p$nage  = p$age_re - p$age_min   #periods before retirement
p$ntr   = p$age_max - p$age_re + 1  #periods after retirement
p$PeriodsToSolve  = p$nage + p$ntr

p$rho     = 2                #(* Coefficient of Relative Risk Aversion *)
p$R       = 1.03             #(* Gross interest rate *)
p$beta    = 0.96             #(* Discount factor *)
p$nP      = 39                #(* Permanent shock Number of points in the discrete approximation to lognormal dist *)
p$nT      = 19                #(* Transitory shock Number of points in the discrete approximation to lognormal dist *)
p$sig2P   = 0.0218           #(* Permanent shock variance of lognormal distribution *)
p$sig2T   = 0.0910           #(* Transitory shock variance of lognormal distribution *)
#p$pUnemp  = 0 #0.5/100      #Probability of unemployment (when unemployed inc level is zero) 
p$aMin  = 1e-6            #(* Lower bound for GothicAVec *)
p$aMax  = 30            #(* Maximum point in GothicAVec *)
p$aHuge = 9000 
p$n = 50                    #(* Number of points in GothicAVec *)
p$NumOfPeople  = 9750     # Number of people to simulate
p$NumOfPeriodsToSimulate = p$PeriodsToSolve   #Length of life in simulation (simulate until age 60)

start_time = proc.time()[3]  
model  <- comp.solveModel(p)
cat(paste('\ntotal seconds to compute Cons rule: ' , proc.time()[3] -  start_time ))

start_time = proc.time()[3]  
moments <- comp.moments(p, model) 
cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))

if(p$age_min==30){
  save(model, moments,p,file='even_can.dat')  
}else{
  save(model, moments,p,file='odd_can.dat') 
} 
