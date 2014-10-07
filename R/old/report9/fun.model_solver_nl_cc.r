source('inc.modelsolver_nl_cc.r')

comp.ngpm <- function(im, p, model){
	source('inc.modelsolver_nl_cc.r')
	
	Mea <- array(0, dim=c(p$nbin, p$ngpa+1))
	Cea <- Mea  
	for ( e in 1:p$nbin){
		Vap=0
		for ( i in 1:p$neps ){
			for ( j in (model$mineta[p$l+1,e] + 1) : (model$maxeta[p$l+1,e] - 1) ){
				lnextm = (p$l * model$mgrid[p$l,im] + model$ygrid[p$l+1,j,i]) / (p$l+1)
				im2 = FindLinProb1(lnextm, model$mgrid[p$l+1,])
				Vap  = Vap + 1/p$neps * model$etaprob[p$l+1,e,j] * 
					( GothVP(p, model$agrid, model$ygrid[p$l+1,j,i], model$C[p$l+1,im2[1],  j,], model$M[p$l+1,im2[1],  j,]) * im2[2] + 
					  GothVP(p, model$agrid, model$ygrid[p$l+1,j,i], model$C[p$l+1,im2[1]+1,j,], model$M[p$l+1,im2[1]+1,j,]) * im2[3] )  # Gothic Va prime	
			}
		}			
		ChiVec = (p$R * p$beta * Vap)^(-1/p$rho) # inverse Euler equation
		MuVec  = model$agrid + ChiVec

		Mea[e,]  = c(0, MuVec)    # Matrix of interpolation data
		Cea[e,]  = c(0,ChiVec)          # Matrix of interpolation data
	}
	return(rbind(Mea,Cea))
}

comp.solveModel <- function(p) {
	eta = comp.eta.prob(p)
	save_eta_name <- paste('eta',p$age_min,'.dat',sep='')	
		
	with( eta, save(xeta, etaprob, mineta, maxeta, etacontot, etauntot, file=save_eta_name) )
	load(save_eta_name)

	#eps: grid and distri = 1/neps
	xeps = comp.eps(p, p$neps)	

	## Profile of income level 
	Incomec = c(1.0304073e+001,1.0343448e+001,1.0382479e+001,1.0421080e+001,1.0459171e+001,1.0496674e+001,1.0533513e+001,1.0569617e+001,1.0604917e+001,1.0639347e+001,1.0672845e+001,1.0705351e+001,1.0736808e+001,1.0767163e+001,1.0796365e+001,1.0824368e+001,1.0851126e+001,1.0876600e+001,1.0900750e+001,1.0923542e+001,1.0944944e+001,1.0964928e+001,1.0983466e+001,1.1000538e+001,1.1016122e+001,1.1030203e+001,1.1042768e+001,1.1053805e+001,1.1063308e+001,1.1071273e+001,1.1077697e+001,1.1082585e+001,1.1085939e+001,1.1087769e+001,1.1088086e+001,1.1086903e+001,1.1084239e+001,1.1080114e+001,1.1074551e+001,1.1067577e+001,1.1059222e+001,1.0369837e+001,1.0369227e+001,1.0368617e+001,1.0368007e+001,1.0367397e+001,1.0366786e+001,1.0366176e+001,1.0365566e+001,1.0364956e+001,1.0364345e+001,1.0363735e+001,1.0363125e+001,1.0362515e+001,1.0361904e+001,1.0361294e+001,1.0360684e+001,1.0360074e+001,1.0359463e+001,1.0358853e+001,1.0358243e+001,1.0357633e+001,1.0357023e+001,1.0356412e+001,1.0355802e+001,1.0355192e+001)
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
	for (it in 1:p$nage){
		vm <- seq(1/(2*p$ngpm), (2*p$ngpm-1)/(2*p$ngpm), l=p$ngpm) #median of each bin
    mgrid[it,] = quantile( yavsim[,it], vm, names=F, na.rm = T )  
	}

	pgrid= mgrid[p$nage,]*0.5 #pension

	agrid = array(0, p$ngpa)
	for( ia in 2:(p$ngpa-1) ){
	  agrid[ia] = agrid[ia-1] + exp(p$pexpgrid*(ia-1))
  }
  ltemp1 = (p$amax-agrid[1])/(agrid[p$ngpa-1]-agrid[1])
  agrid = agrid[1] + ltemp1*(agrid-agrid[1])
  agrid[p$ngpa] = p$aHuge

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
		M[p$nage,im,,]  = array( rep(c(0, MuVec),each=p$nbin),dim=c(p$nbin,p$ngpa+1) )    # Matrix of interpolation data
		C[p$nage,im,,]  = array( rep(c(0,ChiVec),each=p$nbin),dim=c(p$nbin,p$ngpa+1) )    # Matrix of interpolation data
	}

	model= list(	
		mineta= mineta,
		maxeta= maxeta,
		etaprob=etaprob,
		xeps  = xeps,
	  ygrid = ygrid,
	  mgrid = mgrid,
	  pgrid = pgrid,
    agrid = agrid,
	  Mret  = Mret,
	  Cret  = Cret,  
    M    = M,
    C    = C) 

	# lopp over nage-1 ~ 1

	require(snow)  
	cl <- makeCluster(type='MPI')
	chainN = length(cl) 
	cat('Number of Chains: ',chainN,'\n')
	ai = 1:p$ngpm
	
	for ( l in (p$nage-1):1 ){
		p$l <- l

		start_time = proc.time()[3] 
		vals <- parLapply(cl,ai,comp.ngpm,p,model)
		cat(paste('\ntotal seconds for age ', l, ' is', proc.time()[3] -  start_time ))
		
		for(im in ai){
			model$M[l,im,,] <- vals[[im]][ai,]
			model$C[l,im,,] <- vals[[im]][ai+p$ngpm,]
		}	
	}	

	stopCluster(cl)

  return(model)
}

comp.moments<- function(p, model) {
  # This file runs simulation
	ygrid = model$ygrid
  mgrid = model$mgrid
  pgrid = model$pgrid
  agrid = model$agrid
  Mret  = model$Mret
  Cret  = model$Cret 
  M = model$M
  C = model$C
  xeps  = model$xeps

	save_eta_name <- paste('eta',p$age_min,'.dat',sep='')	
	load(save_eta_name)
  set.seed(77)

	# declare variable 
	etasimI <- array( 0,dim=c(p$nsim,p$nage) )
  epssimI <- etasimI
  etasim  <- etasimI
  epssim  <- etasimI
  yavsim	<- etasimI
  msimI   <- etasimI
  msim    <- etasimI
  ysim    <- array( 0,dim=c(p$nsim,p$nage+p$Tret) )
  asim    <- ysim
  xsim    <- ysim
  csim    <- ysim
  #simulate to get distribution of average pre-tax incomes
  etasimI[,1] = sample(1:p$nbin,p$nsim,replace=T)
  epssimI[,1] = sample(1:p$neps,p$nsim,replace=T)

  etasim[,1] = xeta[ cbind(1,etasimI[,1]) ]
  epssim[,1] = xeps[ cbind(1,epssimI[,1]) ]
  ysim[,1]   = ygrid[cbind(1,etasimI[,1],epssimI[,1]) ]
  yavsim[,1]= ysim[,1] 

	for (i in 1:p$nsim){	
		im2 = FindLinProb1(yavsim[i,1],mgrid[1,])
    msimI[i,1] = sample( c(im2[1],im2[1]+1), 1, prob=im2[2:3] )
	}
	msim[,1] = mgrid[ cbind(1,msimI[,1]) ]

  #initial assets
	InitialWYRatio     = c(.17, .5, .83) 
	itemp = sample(1:3,p$nsim,replace=T)
	asim[,1] = InitialWYRatio[itemp] * ysim[,1]

  xsim[,1] = asim[,1] + ysim[,1]
  for (i in 1:p$nsim){
    csim[i,1] = approx( M[1,msimI[i,1],etasimI[i,1],], 
                        C[1,msimI[i,1],etasimI[i,1],], 
                        xsim[i,1] )$y
  }
    
	#loop over life cycle
	for (it in 2:p$nage) {  
    epssimI[,it] = sample(1:p$neps,p$nsim,replace=T)
    epssim[,it]  = xeps[ cbind(it,epssimI[,it]) ]

    etasimI[,it] = sapply(1:p$nsim, function(x)
    	sample(1:p$nbin, 1, prob=etaprob[it, etasimI[x,it-1],]) )
    etasim[,it]  = xeta[ cbind(it,etasimI[,it]) ]
    ysim[,it]    = ygrid[cbind(it,etasimI[,it],epssimI[,it]) ]
    yavsim[,it]  = ( (it-1)*yavsim[,it-1] + ysim[,it] )/it

    for (i in 1:p$nsim){
      im2 = FindLinProb1(yavsim[i,it], mgrid[it,])
      msimI[i,it] = sample( c(im2[1],im2[1]+1), 1, prob=im2[2:3] )
      msim[i,it] = mgrid[it, msimI[i,it]]
      asim[i,it] = p$R*( xsim[i,it-1] - csim[i,it-1] )
      xsim[i,it] = asim[i,it] + ysim[i,it]
      csim[i,it] = approx( M[it, msimI[i,it], etasimI[i,it],], 
                           C[it, msimI[i,it], etasimI[i,it],], 
                           xsim[i,it] )$y     
    }  
  }

##RETIRE####
  for (it in 1:p$Tret) {  
  	ysim[,p$nage+it]  = pgrid[ msimI[,p$nage] ]
  	asim[,p$nage+it]  = p$R * ( xsim[,p$nage+it-1] - csim[,p$nage+it-1] )
  	xsim[,p$nage+it]  = asim[,p$nage+it] + ysim[,p$nage+it]

    for (i in 1:p$nsim){
      csim[i,p$nage+it] = approx( Mret[it, msimI[i,p$nage],], 
                                  Cret[it, msimI[i,p$nage],], 
                                  xsim[i,p$nage+it] )$y
    }  
  }


	moments= list(
  etasimI   = etasimI,
  epssimI   = epssimI,
	etasim    = etasim,
	epssim    = epssim, 
  msimI     = msimI,
  msim      = msim,
  yavsim    = yavsim,
	ysim      = ysim ,
	asim      = asim ,
	xsim      = xsim , 
	csim      = csim)   

	savename <- paste('cohort',p$age_min,'.dat',sep='')
	save(model,moments,p,file=savename)  

  return(moments)
}