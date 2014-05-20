rm(list = ls())
require(EQL)
require(data.table)
require(ggplot2)

setwd('~/git/abbg/R')
source('fun.model_solver_nl_sim.r')

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
#age
p$age_min = 31
p$age_re  = p$age_min+36 #first period income drop 
p$age_max = p$age_min+50
p$nage  = (p$age_re - p$age_min)/2   #periods before retirement
p$ntr   = (p$age_max - p$age_re)/2 + 1  #periods after retirement

#income node
p$nbin  = 100
p$neps  = 23

#asset
p$aMin  = 1e-6            #(* Lower bound for GothicAVec *)
p$aMax  = 30            #(* Maximum point in GothicAVec *)
p$aHuge = 9000 
p$n = 50  

#utility
p$R       = 1.06             #(* Gross interest rate *)
p$beta    = 0.93             #(* Discount factor *)
p$rho     = 2                #(* Coefficient of Relative Risk Aversion *)

#sim
p$nsim  = 9750     # Number of people to simulate
p$N     = 999999

set.seed(77)
start_time = proc.time()[3]  
model  <- comp.solveModel(p)
cat(paste('\ntotal seconds to compute Cons rule: ' , proc.time()[3] -  start_time ))

start_time = proc.time()[3]  
moments <- comp.moments(p, model) 
cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))

if(p$age_min==30){
  save(model, moments,p,file='even.dat')  
}else{
  save(model, moments,p,file='odd.dat') 
} 


#source('fun.sim.data.r')
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
#ddply(sim, ~age,summarise,meanc=mean(lconage),meany=mean(Y))

#persis <- sim.persis(p,sim)
