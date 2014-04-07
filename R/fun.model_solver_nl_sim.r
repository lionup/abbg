source('inc.modelsolver_nl.r')

comp.solveModel <- function(p) {
	res <- with(p,{
		
		#eta: xeta, etaprob, etacontot, etauntot
		#eta <- comp.eta.prob(p,999999)
		#with( eta, save(ieta, xeta, etaprob, etacontot, etauntot, file='eta.dat') )
		load('eta.dat')
		
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
	  load('eta.dat')

		# declare variable 
		epsList  = matrix(0, nrow=nage, ncol=nsim)
		etaList  = epsList
		ctList   = epsList
		stList   = epsList
		mtList   = epsList
		ytList   = epsList
		randeta  = epsList
		enode    = epsList
		visite   = array(0, dim=c(nage, nbin))
	 	
	  # Construct grid for trans draw
	  epsdraws = comp.eps(p, nsim)

		# Construct uniform grid for eta
		etasim <- (1:nsim) / (1+nsim)

		# Initial wy ratio and prob
		InitialWYRatio     = c(.17, .5, .83) 
		InitialWYRatioProb = c(.33333, .33333, .333334)   

		# Construct wtIndicator (list of indicators for initial wealth)
		stIndicator = rep(3, nsim)
		randa <- runif(nsim)
		stIndicator[randa < InitialWYRatioProb[1]] = 1
		stIndicator[randa >= InitialWYRatioProb[1] & randa < (InitialWYRatioProb[1]+InitialWYRatioProb[2])] = 2   

		#loop over life cycle
		for (t in 1:nage) {  
			# Sample randomly from eps grid to get current eps
			epsList[t,] = sample(epsdraws[t,]) 

			# generate random draw on unit interval for current eta
			randeta[t,] = sample(etasim)

			#loop for different individuals
			for (i in 1:nsim){
				# find persistent income draw
		  	if(t==1){  
          enode[1,i] <- which( randeta[1,i] <= etauntot )[1]
          visite[t,enode[1,i]] = visite[t,enode[1,i]] + 1
        }else{
	        enode[t,i] <- which( randeta[t,i] <= etacontot[t,enode[t-1,i],] )[1]
	        visite[t,enode[t,i]] = visite[t,enode[t,i]] + 1
        }

		    etaList[t,i] <- xeta[t,enode[t,i]]

		    # get income
				ytList[t,i] = exp(etaList[t,i] + epsList[t,i])

				if(t==1){
			  	# Construct initial asset
					stList[1,i] = InitialWYRatio[ stIndicator[i] ] * exp(etaList[1,i])
					mtList[1,i] = stList[1,i] + ytList[1,i]      # mtList      : list of normalized m (cash on hand)
				}else{		
	        stList[t,i] = R*( mtList[t-1,i]-ctList[t-1,i] )
	        mtList[t,i] = stList[t,i] + ytList[t,i]
				}
        ctList[t,i] = approx( M[t,enode[t,i],], C[t,enode[t,i],], mtList[t,i] )$y

		  } 
		} 


		model= list(
		epsList    = epsList,
		etaList    = etaList, 
		ctList     = ctList ,
		stList     = stList ,
		mtList     = mtList , 
		ytList     = ytList  , 
		enode      = enode  ,
		visite     = visite )   
	}) 

  return(res)
}