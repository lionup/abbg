source('inc.modelsolver_nl.r')

comp.solveModel <- function(p, eta) {
	res <- with(p,{
		
		#eta: xeta, etaprob, etacontot, etauntot
		#eta <- comp.eta(p)
		#with( eta, save(xeta, etaprob, etacontot, etauntot, file='eta.dat') )
		load('eta.dat')
		ieta = exp(xeta)
		
		#eps: grid and distri 
		epsprob = rep(1/neps,neps)
		xeps = comp.eps(p)	
		ieps = exp(xeps) 

		#not being able to borrow more than can repay for sure 
		#this is the lower bound of asset at the end of the period
 		mininc <- rep(0, nage) # in the last period, end asset is 0
		for(t in (nage-1):1){
	    mininc[t] = ( mininc[t+1] + ieta[t+1,1] * ieps[t+1,1] ) / R 
		}	

		#asset grid
    AlphaVec=exp(exp(exp(seq(log(log(log(aMin+1)+1)+1),log(log(log(aMax+1)+1)+1),l=n))-1)-1)-1  # set up triple exponential grid
		AlphaVec = c(AlphaVec,aHuge)

		# Construct array of interpolation data
		C = array( rep(0:(n+1), each=nage*nbin), dim=c(nage, nbin, n+2) )
		M = array( rep(0:(n+1), each=nage*nbin), dim=c(nage, nbin, n+2) )

		# lopp over T-1 ~ 1
		for ( l in (nage-1):1 ){	
			for ( e in 1:nbin){
				# end asset grid this period
				AlphaVec1 = AlphaVec - mininc[l]

				# eta this period ieta[l,e]
				# get eta next period ieta[l+1,]
				# get the distri of eta next period etaprob[l+1,e,]
				#get eps next period ieps[l+1,]
				#get the distri of eps next period	epsprob
				Vap    = GothVP(p, AlphaVec1, ieps[l+1,], ieta[l+1,], 
					epsprob, etaprob[l+1,e,], C[l+1,,], M[l+1,,] )   # Gothic Va prime
				ChiVec = (R * beta * Vap)^(-1/rho) # inverse Euler equation
	 			MuVec  = AlphaVec1+ChiVec
	  		M[l,e,]  = c(-mininc[l], MuVec)    # Matrix of interpolation data
	  		C[l,e,]  = c(0,ChiVec)          # Matrix of interpolation data
			}
		}

		model= list(
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