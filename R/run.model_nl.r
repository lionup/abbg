rm(list = ls())
require(EQL)
require(data.table)
require(ggplot2)

setwd('~/git/abbg/R')
source('fun.model_solver_nl_sim.r')
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

detach(data)
#age
p$age_min = 31
p$age_re  = p$age_min+36 #income drop almost half
p$age_max = p$age_min+50
p$nage  = (p$age_re - p$age_min)/2   #periods before retirement
p$ntr   = (p$age_max - p$age_re)/2 + 1  #periods after retirement

#income node
p$nbin  = 23
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
p$nsim  = 10000     # Number of people to simulate

set.seed(77)

start_time = proc.time()[3]  
model  <- comp.solveModel(p)
cat(paste('\ntotal seconds to compute Cons rule: ' , proc.time()[3] -  start_time ))

start_time = proc.time()[3]  
moments <- comp.moments(p, model) 
cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))

age = seq(p$age_min, p$age_max, 2)
mm <- data.frame( age=age)

attach(moments)
#median
mm <- cbind( mm, data.frame(stMedian = apply(stList, 1, median)) )
mm <- cbind( mm, data.frame(ctMedian = apply(ctList, 1, median)) )
mm <- cbind( mm, data.frame(mtMedian = apply(mtList, 1, median)) )
mm <- cbind( mm, data.frame(ytMedian = apply(ytList, 1, median)) )

#1st quantile
mm <- cbind( mm, data.frame(st1q = apply(stList, 1, quantile, 0.25)) )
mm <- cbind( mm, data.frame(ct1q = apply(ctList, 1, quantile, 0.25)) )
mm <- cbind( mm, data.frame(mt1q = apply(mtList, 1, quantile, 0.25)) )
mm <- cbind( mm, data.frame(yt1q = apply(ytList, 1, quantile, 0.25)) )

#3st quantil
mm <- cbind( mm, data.frame(st3q = apply(stList, 1, quantile, 0.75)) )
mm <- cbind( mm, data.frame(ct3q = apply(ctList, 1, quantile, 0.75)) )
mm <- cbind( mm, data.frame(mt3q = apply(mtList, 1, quantile, 0.75)) )
mm <- cbind( mm, data.frame(yt3q = apply(ytList, 1, quantile, 0.75)) )

detach(moments)

mmic <- data.frame( age=age, value = mm$yt1q, moment = 'income', quantile='1st' )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ytMedian, 
	moment = 'income', quantile='Median') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$yt3q, 
	moment = 'income', quantile='3rd') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ct1q, 
	moment = 'consumption', quantile='1st') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ctMedian, 
	moment = 'consumption', quantile='Median') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ct3q, 
	moment = 'consumption', quantile='3rd') )

mma <- data.frame( age=age, asset = mm$st1q, quantile='1st') 
mma <- rbind(mma, data.frame( age=age, asset = mm$stMedian, quantile='Median') )
mma <- rbind(mma, data.frame( age=age, asset = mm$st3q, quantile='3rd') )


save(model, moments,mm,file='1e5.dat')