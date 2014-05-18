source('inc.modelsolver_nl.r')

comp.solveModel <- function(p) {
	res <- with(p,{
		
		#eta: xeta, etaprob, etacontot, etauntot
		eta <- comp.eta.prob(p)
		if(age_min==30){
			with( eta, save(ieta, xeta, etaprob, etacontot, etauntot, file='eta_even.dat') )
			load('eta_even.dat')
		}else{
			with( eta, save(ieta, xeta, etaprob, etacontot, etauntot, file='eta_odd.dat') )
			load('eta_odd.dat')
		}	

		#eps: grid and distri 
		epsprob = rep(1/neps,neps)
		xeps = comp.eps(p, neps)	
		ieps = exp(xeps)
		#eps after retirement
		ieps = rbind( ieps,matrix(1,nrow=ntr, ncol=neps) )

		## Profile of income level 
		Incomec = c(1.0304073e+001,1.0343448e+001,1.0382479e+001,1.0421080e+001,1.0459171e+001,1.0496674e+001,1.0533513e+001,1.0569617e+001,1.0604917e+001,1.0639347e+001,1.0672845e+001,1.0705351e+001,1.0736808e+001,1.0767163e+001,1.0796365e+001,1.0824368e+001,1.0851126e+001,1.0876600e+001,1.0900750e+001,1.0923542e+001,1.0944944e+001,1.0964928e+001,1.0983466e+001,1.1000538e+001,1.1016122e+001,1.1030203e+001,1.1042768e+001,1.1053805e+001,1.1063308e+001,1.1071273e+001,1.1077697e+001,1.1082585e+001,1.1085939e+001,1.1087769e+001,1.1088086e+001,1.1086903e+001,1.1084239e+001,1.1080114e+001,1.1074551e+001,1.1067577e+001,1.1059222e+001,1.0369837e+001,1.0369227e+001,1.0368617e+001,1.0368007e+001,1.0367397e+001,1.0366786e+001,1.0366176e+001,1.0365566e+001,1.0364956e+001,1.0364345e+001,1.0363735e+001,1.0363125e+001,1.0362515e+001,1.0361904e+001,1.0361294e+001,1.0360684e+001,1.0360074e+001,1.0359463e+001,1.0358853e+001,1.0358243e+001,1.0357633e+001,1.0357023e+001,1.0356412e+001,1.0355802e+001,1.0355192e+001)
		#Incomec = rep(0,rep=66)
		incpro = data.table(age=25:90,income=Incomec)
		incpro = subset(incpro,age >= age_min & age <= age_max)
		incpro[,IncomeLevelc := exp(income)]
		setkey(incpro, age)
		age_full = seq(age_min, age_max, 2)
		incpro = incpro[J(age_full)]
		incpro[,norminc:=IncomeLevelc/exp(Incomec[6])]
		# Income growth factors for t+1~T 
		#Glist = with(incpro, { IncomeLevelc[-1]/IncomeLevelc[-length(IncomeLevelc)] })
		inc = incpro$norminc

		#not being able to borrow more than can repay for sure 
		#this is the lower bound of asset at the end of the period
 		mininc <- rep(0, nage+ntr) # in the last period, end asset is 0
		for(t in (nage+ntr-1):1){
	    mininc[t] = ( mininc[t+1]+inc[t+1]*ieta[t+1,1]*ieps[t+1,1] ) / R 
		}	

		#asset grid
    AlphaVec=exp(exp(exp(seq(log(log(log(aMin+1)+1)+1),log(log(log(aMax+1)+1)+1),l=n))-1)-1)-1  # set up triple exponential grid
		AlphaVec = c(AlphaVec,aHuge)

		# Construct array of interpolation data
		C = array( rep(0:(n+1), each=(nage+ntr)*nbin), dim=c((nage+ntr), nbin, n+2) )
		M = array( rep(0:(n+1), each=(nage+ntr)*nbin), dim=c((nage+ntr), nbin, n+2) )

		# lopp over T-1 ~ 1
		for ( l in ((nage+ntr)-1):1 ){	
			for ( e in 1:nbin){
				# end asset grid this period
				AlphaVec1 = AlphaVec - mininc[l]

				# eta this period ieta[l,e]
				# get eta next period ieta[l+1,]
				# get the distri of eta next period etaprob[l+1,e,]
				#get eps next period ieps[l+1,]
				#get the distri of eps next period	epsprob
				Vap    = GothVP(p, AlphaVec1, inc[l+1], ieps[l+1,], ieta[l+1,], 
					epsprob, etaprob[l+1,e,], C[l+1,,], M[l+1,,] )   # Gothic Va prime
				ChiVec = (R * beta * Vap)^(-1/rho) # inverse Euler equation
	 			MuVec  = AlphaVec1+ChiVec
	  		M[l,e,]  = c(-mininc[l], MuVec)    # Matrix of interpolation data
	  		C[l,e,]  = c(0,ChiVec)          # Matrix of interpolation data
			}
		}

		model= list(
	  inc  = inc,
    M    = M,
    C    = C)
	}) 

  return(res)
}

comp.moments<- function(p, model) {
	res <- with(p,{  
	  # This file runs simulation
	  inc = model$inc
	  M = model$M
	  C = model$C

		if(age_min==30){ 
			load('eta_even.dat')
		}else{
			load('eta_odd.dat')
		}

		# declare variable 
		epsList  = matrix(0, nrow=(nage+ntr), ncol=nsim)
		etaList  = epsList
		ctList   = epsList
		stList   = epsList
		mtList   = epsList
		ytList   = epsList
		randeta  = epsList
		enode    = epsList
		visite   = array(0, dim=c((nage+ntr), nbin))
	 	
	 	set.seed(77)
	  # Construct grid for trans draw
	  epsdraws = comp.eps(p, nsim)
		#eps after retirement
		epsdraws = rbind( epsdraws,matrix(0,nrow=ntr, ncol=nsim) )

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
		for (t in 1:(nage+ntr)) {  
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
				ytList[t,i] =  inc[t] * exp( etaList[t,i] + epsList[t,i] )

				if(t==1){
			  	# Construct initial asset
					stList[1,i] = InitialWYRatio[ stIndicator[i] ] * inc[1] * exp(etaList[1,i])
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