source('inc.modelsolver_rw.r')

comp.solveModel <- function(p) {
	res <- with(p,{  
		
		#Setting up shock values (discrete approximation to log normal)   
		PermVecProb = rep(1/nP,nP)
		ThetaVecProb = rep(1/nT,nT)
		#if (pUnemp>0) { # If assume unemployment
   	# ThetaVecProb = c( pUnemp,ThetaVecProb*(1-pUnemp) )
 		#}

		PermVec  = DiscreteApproxToMeanOneLogNormal(sig2P,nP)
		ThetaVec = DiscreteApproxToMeanOneLogNormal(sig2T,nT)
		#if (pUnemp>0){ # If assume unemployment
    #	ThetaVec = c( 0,ThetaVec/(1-pUnemp) )
 		#}
    
    #Theta 
		ThetaMat = t( replicate(nage, ThetaVec) )
		ThetaMat = rbind( ThetaMat, matrix(1,nrow=ntr,ncol=length(ThetaVec)) )

		PermMat = t( replicate(nage, PermVec) )
		PermMat = rbind( PermMat, matrix(1,nrow=ntr,ncol=length(PermVec)) )

		## Profile of income level 
		Incomec = c(1.0304073e+001,1.0343448e+001,1.0382479e+001,1.0421080e+001,1.0459171e+001,1.0496674e+001,1.0533513e+001,1.0569617e+001,1.0604917e+001,1.0639347e+001,1.0672845e+001,1.0705351e+001,1.0736808e+001,1.0767163e+001,1.0796365e+001,1.0824368e+001,1.0851126e+001,1.0876600e+001,1.0900750e+001,1.0923542e+001,1.0944944e+001,1.0964928e+001,1.0983466e+001,1.1000538e+001,1.1016122e+001,1.1030203e+001,1.1042768e+001,1.1053805e+001,1.1063308e+001,1.1071273e+001,1.1077697e+001,1.1082585e+001,1.1085939e+001,1.1087769e+001,1.1088086e+001,1.1086903e+001,1.1084239e+001,1.1080114e+001,1.1074551e+001,1.1067577e+001,1.1059222e+001,1.0369837e+001,1.0369227e+001,1.0368617e+001,1.0368007e+001,1.0367397e+001,1.0366786e+001,1.0366176e+001,1.0365566e+001,1.0364956e+001,1.0364345e+001,1.0363735e+001,1.0363125e+001,1.0362515e+001,1.0361904e+001,1.0361294e+001,1.0360684e+001,1.0360074e+001,1.0359463e+001,1.0358853e+001,1.0358243e+001,1.0357633e+001,1.0357023e+001,1.0356412e+001,1.0355802e+001,1.0355192e+001)
		incpro = data.table(age=25:90,income=Incomec)
		incpro = subset(incpro,age >= age_min & age <= age_max)
		incpro[,IncomeLevelc := exp(income)]
		setkey(incpro, age)
		age_full = seq(age_min, age_max, 2)
		incpro = incpro[J(age_full)]
		incpro[,norminc:=IncomeLevelc/exp(Incomec[6])]
		inc = incpro$norminc

		# Probability of being alive after retirement 
		# (1st element is the prob of being alive until age 66)
		#ProbOfAlive = c(9.8438596e-01,9.8438596e-01,9.8438596e-01,9.8438596e-01,9.8438596e-01,9.7567062e-01,9.7567062e-01,9.7567062e-01,9.7567062e-01,9.7567062e-01,9.6207901e-01,9.6207901e-01,9.6207901e-01,9.6207901e-01,9.6207901e-01,9.3721595e-01,9.3721595e-01,9.3721595e-01,9.3721595e-01,9.3721595e-01,6.3095734e-01,6.3095734e-01,6.3095734e-01,6.3095734e-01,6.3095734e-01)
		#ProbOfAlive = c(rep(1,PeriodsToSolve-length(ProbOfAlive)),ProbOfAlive)

		# Corrected beta: Exp(dZ(t)*theta)
	  #Betacorr = c(1.0649141e+00,1.0579968e+00,1.0514217e+00,1.0451790e+00,1.0392591e+00,1.0336529e+00,1.0283519e+00,1.0233477e+00,1.0186323e+00,1.0141979e+00,1.0100373e+00,1.0061433e+00,1.0025092e+00,9.9912824e-01,9.9599427e-01,9.9310116e-01,9.9044306e-01,9.8801430e-01,9.8580946e-01,9.8382325e-01,9.8205060e-01,9.8048658e-01,9.7912645e-01,9.7796558e-01,9.7699952e-01,9.7622393e-01,9.7563459e-01,9.7522741e-01,9.7499842e-01,9.7494373e-01,9.7505955e-01,9.7534220e-01,9.7578805e-01,9.7639357e-01,9.7715529e-01,9.7806981e-01,9.7913379e-01,9.8034393e-01,9.8169700e-01,8.2872135e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01,9.9021110e-01)

		#asset grid
    AlphaVec=exp(exp(exp(seq(log(log(log(aMin+1)+1)+1),log(log(log(aMax+1)+1)+1),l=n))-1)-1)-1  # set up triple exponential grid
		AlphaVec = c(AlphaVec,aHuge)

		# Construct matrix of interpolation data
		#	Period T, min of cash is 0, so min of cons is 0
 		C = matrix(0:(n+1))
		M = matrix(0:(n+1))
		PDVmwn = 0 # present disc value of next period min wage at T

		# lopp over T-1 ~ 1
		for (l in 1:PeriodsToSolve){	
			# calculate ct from each grid point in AlphaVec
			#P      =  R * GList[65-l+1]^(-rho) * beta * ProbOfAlive[length(ProbOfAlive)-l+1] * Betacorr[length(Betacorr)-l+1]
			G      = GList[65-l+1]
			ThetaVec = ThetaMat[65-l+1,]
			PermVec  = PermMat[65-l+1,]

			PDVmwn = ( PDVmwn + inc[l]*min(ThetaVec)*min(PermVec) )/R
			AlphaVec1=AlphaVec-PDVmwn
			Vap    = GothVP(p, AlphaVec1, G, ThetaVec, PermVec, ThetaVecProb, PermVecProb, C, M )   # Gothic Va prime
			ChiVec = (R * beta * G^(-rho) * Vap)^(-1/rho) # inverse Euler equation
 			MuVec  = AlphaVec1+ChiVec
  		M      = cbind(M, c(-PDVmwn, MuVec))                  # Matrix of interpolation data
  		C      = cbind(C, c(0,ChiVec))                 # Matrix of interpolation data
		}

		model= list(
		GList= GList,
    M    = M,
    C    = C)
	}) 

  return(res)
}

comp.moments<- function(p, model) {
	res <- with(p,{  

	  # This file runs simulation
	  GList= model$GList
	  M = model$M
	  C = model$C

	  #declare variable 
	  ThetaList= matrix(0, nrow=NumOfPeriodsToSimulate, ncol=NumOfPeople)
	  PermList = ThetaList 
		ctList   = ThetaList
		stList   = ThetaList
		mtList   = ThetaList
		Perm     = ThetaList
		Income   = ThetaList

		# Construct income shock draw lists
		# Tran shock draw list
		ThetaDraws = DiscreteApproxToMeanOneLogNormal(sigT,NumOfPeople)
		#if (pUnemp>0){ # If assume unemployment
		#	ThetaDraws = DiscreteApproxToMeanOneLogNormal(sigT,NumOfPeople-pUnemp*NumOfPeople)
		#	ThetaDraws = ThetaDraws/(1 - pUnemp)
		#	ThetaDraws = append(rep(0,pUnemp*NumOfPeople),ThetaDraws)
		#} 

		# Perm shock draw list
		PermShockDraws = DiscreteApproxToMeanOneLogNormal(sigP,NumOfPeople)

		# Construct income shock lists
		for (i in 1:NumOfPeriodsToSimulate){
		  ThetaList[i,] = sample(ThetaDraws)     # List of Theta (tran shock) 
		  PermList[i,]  = sample(PermShockDraws) # List of perm shock 
		}

		if (NumOfPeriodsToSimulate > 41){
			PermList[42:NumOfPeriodsToSimulate,] = matrix(1,nrow=NumOfPeriodsToSimulate-41,ncol=length(PermShockDraws)) 
			ThetaList[42:NumOfPeriodsToSimulate,]= matrix(1,nrow=NumOfPeriodsToSimulate-41,ncol=length(ThetaDraws)) 
		}

		# Construct wtIndicator (list of indicators for initial wealth)
		InitialWYRatio     = c(.17, .5, .83)              # Initial wy ratio (from the program on the paper (p.13))
		InitialWYRatioProb = c(.33333, .33333, .333334)   # Prob associated with initial wy ratio 

		stIndicator = rep(0, NumOfPeople)
		for (i in 1:NumOfPeople){
	    r = runif(1)
	    if (r < InitialWYRatioProb[1]){
	      stIndicator[i] = 1
	    }else if (r < InitialWYRatioProb[1]+InitialWYRatioProb[2]){
	      stIndicator[i] = 2
	    }else{
	      stIndicator[i] = 3
	    }
		}
		
		# First period
		stList[1,] = InitialWYRatio[ stIndicator ] # stList      : list of normalized s (savings at the beginning of age)       
		mtList[1,] = stList[1,] + ThetaList[1,]      # mtList      : list of normalized m (cash on hand)

		# Construct itIndicator (list of indicators for initial income)
		InitialIncome     = c(10000, 20000, 40000)              # Initial wy ratio (from the program on the paper (p.13))
		InitialIncomeProb = c(.3, .4, .3)   # Prob associated with initial wy ratio 

		itIndicator = rep(0, NumOfPeople)
		for (i in 1:NumOfPeople){
	    r = runif(1)
	    if (r < InitialIncomeProb[1]){
	      itIndicator[i] = 1
	    }else if (r < InitialIncomeProb[1]+InitialIncomeProb[2]){
	      itIndicator[i] = 2
	    }else{
	      itIndicator[i] = 3
	    }
		}
		
		Perm[1,] = InitialIncome[ itIndicator ] * PermList[1,]
		Income[1,]= Perm[1,] * ThetaList[1,]

		len = ncol(M)
		ctList[1,] =  approx( M[,len], C[,len], mtList[1,] )$y
		#ctList[1, mtList[1,] == 0] =0                     

		# Continue simulstion
		for (t in 2:NumOfPeriodsToSimulate){
	    for (j in 1:NumOfPeople){           
        stList[t,j] = ( R/GList[t-1]/PermList[t,j] )*( mtList[t-1,j]-ctList[t-1,j] )
        mtList[t,j] = stList[t,j] + ThetaList[t,j]
        ctList[t,j] = approx( M[,len-t+1], C[,len-t+1], mtList[t,j] )$y
        #if (mtList[t,j] == 0) ctList[t,j] =0                       # ctList      : list of normalized consumption 
				Perm[t,j]   = Perm[t-1,j] * GList[t-1] * PermList[t,j] # List of perm income
				Income[t,j] = Perm[t,j] * ThetaList[t,j]
			}
		}

		stMeanList   = rowMeans(stList)         # stLevelMedianList: list of median of savings level 
		stMedianList = apply(stList, 1, median)
		ctMeanList   = rowMeans(ctList) 
		ctMedianList = apply(ctList, 1, median)
		Ct           = ctList * Perm
		CtMean       = rowMeans(Ct)
		IncomeMean   = rowMeans(Income)
		
		model= list(
		stList        = stList,
		stMeanList    = stMeanList,
		stMedianList  = stMedianList, 
		ctList        = ctList, 
		ctMeanList    = ctMeanList,
		ctMedianList  = ctMedianList,
		Perm          = Perm,
		Income        = Income, 
		IncomeMean    = IncomeMean, 
		Ct            = Ct, 
		CtMean        = CtMean)
	}) 

  return(res)
}