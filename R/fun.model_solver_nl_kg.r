source('inc.modelsolver_nl_kg.r')

comp.solveModel <- function(p) {
	#eta = comp.eta.prob(p)
	save_eta_name <- paste('eta',p$age_min,'.dat',sep='')	
		
	#with( eta, save(xeta, etaprob, mineta, maxeta, etacontot, etauntot, file=save_eta_name) )
	load(save_eta_name)

	#eps: grid and distri 
	epsprob = rep(1/p$neps,p$neps)
	xeps = comp.eps(p, p$neps)	

	## Profile of income level 
	Incomec = c(1.0304073e+001,1.0343448e+001,1.0382479e+001,1.0421080e+001,1.0459171e+001,1.0496674e+001,1.0533513e+001,1.0569617e+001,1.0604917e+001,1.0639347e+001,1.0672845e+001,1.0705351e+001,1.0736808e+001,1.0767163e+001,1.0796365e+001,1.0824368e+001,1.0851126e+001,1.0876600e+001,1.0900750e+001,1.0923542e+001,1.0944944e+001,1.0964928e+001,1.0983466e+001,1.1000538e+001,1.1016122e+001,1.1030203e+001,1.1042768e+001,1.1053805e+001,1.1063308e+001,1.1071273e+001,1.1077697e+001,1.1082585e+001,1.1085939e+001,1.1087769e+001,1.1088086e+001,1.1086903e+001,1.1084239e+001,1.1080114e+001,1.1074551e+001,1.1067577e+001,1.1059222e+001,1.0369837e+001,1.0369227e+001,1.0368617e+001,1.0368007e+001,1.0367397e+001,1.0366786e+001,1.0366176e+001,1.0365566e+001,1.0364956e+001,1.0364345e+001,1.0363735e+001,1.0363125e+001,1.0362515e+001,1.0361904e+001,1.0361294e+001,1.0360684e+001,1.0360074e+001,1.0359463e+001,1.0358853e+001,1.0358243e+001,1.0357633e+001,1.0357023e+001,1.0356412e+001,1.0355802e+001,1.0355192e+001)
	#Incomec = rep(0,rep=66)
	incpro = data.table(age=25:90,income=Incomec)
	setkey(incpro, age)
	age_full = seq(p$age_min, p$age_re-2, 2)
	incpro = incpro[J(age_full)]
	kappa = incpro$income

	ygrid = array(0,dim=c(p$nage,p$nbin,p$neps))
	for(it in 1:p$nage){
	  for(iz in 1:p$nbin){
    	for(ie in 1:p$neps){
    		ygrid[it,iz,ie] = exp(kappa[it] + xeta[it,iz] + xeps[it,ie])
    	}
    }
  } 	

  zsimI <- array( 0,dim=c(p$nsim,p$nage) )
  esimI <- zsimI
  ysim  <- zsimI
  yavsim<- zsimI
  #simulate to get distribution of average pre-tax incomes
  zsimI[,1] = sample(1:p$nbin,p$nsim,replace=T)
  esimI[,1] = sample(1:p$neps,p$nsim,replace=T)

	for (i in 1:p$nsim){
		ysim[i,1] = ygrid[ 1,zsimI[i,1],esimI[i,1] ]
		yavsim[i,1]= ysim[i,1] 
	  for (it in 2:p$nage){
	     zsimI[i,it] = sample( 1:p$nbin,1,prob=etaprob[it,zsimI[i,it-1],] )
	     esimI[i,it] = sample( 1:p$neps,1)
	     ysim[i,it]  = ygrid[ it,zsimI[i,it],esimI[i,it] ]
	     yavsim[i,it]= ( (it-1)*yavsim[i,it-1] + ysim[i,it] )/it
		}	 
	}
	        
	mgrid = array(0,dim=c(p$nage,p$ngpm))        
	#equally spaced between 5th perc and median, and median and 95th perc
	for (it in 1:p$nage){
    lysim = quantile( yavsim[,it], c(0.05,0.5,0.95), names=F, na.rm = T )  
    midm = (1+p$ngpm)/2
    mgrid[it,1]    = lysim[1] 
    mgrid[it,midm] = lysim[2] 
    mgrid[it,p$ngpm] = lysim[3] 

    lwidth1 = (mgrid[it,midm]  - mgrid[it,1])   /(midm-1)
    lwidth2 = (mgrid[it,p$ngpm]- mgrid[it,midm])/(midm-1)
    
    for(im in 2:(midm-1)){
      mgrid[it,im] = mgrid[it,im-1] +lwidth1
    }
    for( im in (midm+1):(p$ngpm-1) ) {
      mgrid[it,im] = mgrid[it,im-1] +lwidth2
    }
	}

	pgrid= mgrid[p$nage,]*0.5 #pension

	agrid = array(0, p$ngpa)
	for(ia in 2:p$ngpa){
	  agrid[ia] = agrid[ia-1] + exp(p$pexpgrid*(ia-1))
  }
  ltemp1 = (p$amax-agrid[1])/(agrid[p$ngpa]-agrid[1])
  agrid = agrid[1] + ltemp1*(agrid-agrid[1])

  # create last period C and M after retirement
 	Cret = array( rep(seq(0,1e7,l=p$ngpa+1), each=p$Tret*p$ngpm), 
 		dim=c(p$Tret,p$ngpm,p$ngpa+1) )
	Mret = Cret

	# lopp over retirement
	for ( l in (p$Tret-1):1 ){	
		for ( im in 1:p$ngpm){
			Vap    = GothVP(p, agrid, pgrid[im], Cret[l+1,im,], Mret[l+1,im,])   # Gothic Va prime
			ChiVec = (p$R * p$beta * Vap)^(-1/p$rho) # inverse Euler equation
 			MuVec  = agrid + ChiVec
  		Mret[l,im,]  = c(0, MuVec)    # Matrix of interpolation data
  		Cret[l,im,]  = c(0,ChiVec)          # Matrix of interpolation data
		}
	}

	# create C and M before retirement
 	C = array( 0, dim=c(p$nage,p$ngpm,p$nbin,p$ngpa+1) )
	M = C

	#last period of working
	for ( im in 1:p$ngpm){
		Vap    = GothVP(p, agrid, pgrid[im], Cret[1,im,], Mret[1,im,])   # Gothic Va prime
		ChiVec = (p$R * p$beta * Vap)^(-1/p$rho) # inverse Euler equation
		MuVec  = agrid + ChiVec
		M[p$nage,im,,]  = rep(c(0, MuVec),each=p$nbin)    # Matrix of interpolation data
		C[p$nage,im,,]  = rep(c(0,ChiVec),each=p$nbin)    # Matrix of interpolation data
	}

	# lopp over nage-1 ~ 1
	for ( l in (p$nage-1):1 ){
		for ( im in 1:p$ngpm){	
			for ( e in 1:p$nbin){
				Vap=0
    		for ( i in 1:p$neps ){
    			for ( j in (mineta[l+1,e]+1):(maxeta[l+1,e]-1) ){
    				lnextm = (l*mgrid[l,im] + ygrid[l+1,j,i])/(l+1)
    				im2 = FindLinProb1(lnextm,mgrid[l+1,])
						Vap  = Vap + epsprob[i] * etaprob[l+1,e,j] * 
							( GothVP(p, agrid, ygrid[l+1,j,i], C[l+1,im2[1],  j,], M[l+1,im2[1],  j,]) * im2[2] + 
							  GothVP(p, agrid, ygrid[l+1,j,i], C[l+1,im2[1]+1,j,], M[l+1,im2[1]+1,j,]) * im2[3] )  # Gothic Va prime	
					}
				}			
				ChiVec = (p$R * p$beta * Vap)^(-1/p$rho) # inverse Euler equation
	 			MuVec  = agrid+ChiVec
	  		M[l,im,e,]  = c(0, MuVec)    # Matrix of interpolation data
	  		C[l,im,e,]  = c(0,ChiVec)          # Matrix of interpolation data
			}
		}
	}		

	res= list(
	  ygrid = ygrid,
	  mgrid = mgrid,
	  pgrid = pgrid,
	  Mret  = Mret,
	  Cret  = Cret,  
    M    = M,
    C    = C) 

  return(res)
}

comp.moments<- function(p, model) {
  # This file runs simulation
	ygrid = model$ygrid
  mgrid = model$mgrid
  pgrid = model$pgrid
  Mret  = model$Mret
  Cret  = model$Cret 
  M = model$M
  C = model$C

	save_eta_name <- paste('eta',p$age_min,'.dat',sep='')	
	load(save_eta_name)

	# declare variable 
	epsList  = matrix(0, nrow=p$nage, ncol=p$nsim)
	etaList  = epsList
	ctList   = epsList
	stList   = epsList
	mtList   = epsList
	ytList   = epsList
	randeta  = epsList
	enode    = epsList
 	
 	set.seed(77)
  # Construct grid for trans draw
  epsdraws = comp.eps(p, p$nsim)
	etasim <- (1:p$nsim) / (1+p$nsim)

	etasimI <- array( 0,dim=c(p$nsim,p$nage) )
  epssimI <- etasimI
  ysim  	<- etasimI
  yavsim	<- etasimI
  #simulate to get distribution of average pre-tax incomes
  etasimI[,1] = sample(1:p$nbin,p$nsim,replace=T)
  epssimI[,1] = sample(1:p$neps,p$nsim,replace=T)

	for (i in 1:p$nsim){
		ysim[i,1] = ygrid[ 1,zsimI[i,1],esimI[i,1] ]
		yavsim[i,1]= ysim[i,1] 
		im2 = FindLinProb1(yavsim[i,1],mgrid[1,])
	}



	# Initial wy ratio and prob
	InitialWYRatio     = c(.17, .5, .83) 
	InitialWYRatioProb = c(.33333, .33333, .333334)   

	# Construct wtIndicator (list of indicators for initial wealth)
	stIndicator = rep(3, nsim)
	randa <- runif(nsim)
	stIndicator[randa < InitialWYRatioProb[1]] = 1
	stIndicator[randa >= InitialWYRatioProb[1] & randa < (InitialWYRatioProb[1]+InitialWYRatioProb[2])] = 2   

	#loop over life cycle
	for (t in 1:p$nage) {  
		# Sample randomly from eps grid to get current eps
		epsList[t,] = sample(epsdraws[t,]) 

		# generate random draw on unit interval for current eta
		randeta[t,] = sample(etasim)

		#loop for different individuals
		for (i in 1:p$nsim){
			# find persistent income draw
	  	if(t==1){  
        enode[1,i] <- which( randeta[t,i] <= etauntot )[1]
      }else{
        enode[t,i] <- which( randeta[t,i] <= etacontot[t,enode[t-1,i],] )[1]
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
	enode      = enode )   

	savename <- paste('cohort',p$age_min,'.dat',sep='')
	save(model,p,file=savename)  

  return(model)
}

comp.income <- function(p){
	res <- with(p,{

		eta <- comp.eta.prob(p)
		save_eta_name <- paste('eta',age_min,'.dat',sep='')	
		with( eta, save(ieta, xeta, etaprob, mineta, maxeta, etacontot, etauntot, file=save_eta_name) )
		load(save_eta_name)

		epsList  = matrix(0, nrow=(nage+ntr), ncol=nsim)
		etaList  = epsList
		randeta  = epsList
		enode    = epsList
		visite   = array(0, dim=c((nage+ntr), nbin))
	 	
	 	set.seed(77)
	  # Construct grid for trans draw
	  epsdraws = comp.eps(p, nsim)
		#eps after retirement
		epsdraws = rbind( epsdraws,matrix(0,nrow=ntr, ncol=nsim) )

		etasim <- (1:nsim) / (1+nsim)

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
		  } 
		} 


		model= list(
		epsList    = epsList,
		etaList    = etaList, 
		enode      = enode  ,
		visite     = visite )   
	}) 

  return(res)

}