source('inc.modelsolver_nl.r')

comp.solveModel <- function(p) {
	res <- with(p,{
		
		#eta: xeta, etaprob, etacontot, etauntot
		eta <- comp.eta.prob(p,999999)
		with( eta, save(xeta, etaprob, etacontot, etauntot, file='eta.dat') )
		load('eta.dat')
		ieta = exp(xeta)
		
		#eps: grid and distri 
		epsprob = rep(1/neps,neps)
		xeps = comp.eps(p, neps)	
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
	  M = model$M
	  C = model$C

		# declare variable 
		epsList= matrix(0, nrow=nage, ncol=nsim)
		ctList   = epsList
		stList   = epsList
		mtList   = epsList

	 	# Construct transitory income shock draw lists
	  # Construct grid for trans draw
	  epsdraws = comp.eps(p, nsim)

	  # Sample randomly from the grid to get the distri for each period
	  for (i in 1:nage){
	  	epsList[i,] = sample(epsdraws[i,])
	  }

		# Construct persistent income shock draw lists
		etaList = comp.eta.sim(p, nsim)

		# get income
		yList = exp(etaList + epsList)

		# Construct initial asset
		# Construct wtIndicator (list of indicators for initial wealth)
		InitialWYRatio     = c(.17, .5, .83)              # Initial wy ratio (from the program on the paper (p.13))

		stIndicator = rep(3, nsim)
		rama <- runif(nsim)
		stIndicator[rama < InitialWYRatioProb[1]] = 1
		stIndicator[rama < InitialWYRatioProb[1]+InitialWYRatioProb[2]] = 2

		# First period
		stList[1,] = InitialWYRatio[ stIndicator ] # stList      : list of normalized s (savings at the beginning of age)       
		mtList[1,] = stList[1,] + yList[1,]      # mtList      : list of normalized m (cash on hand)

		# first period consumption

		# Simulate


	}) 

  return(res)
}