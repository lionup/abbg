
comp.solveModel <- function(p) {
	res <- with(p,{ 
		#####################################################
		#DATA
		#####################################################
		#log experience profile
		kappa = read.table("~/git/abbg/R/old/KV/Input/kappasmooth.txt")$V1 #25~59
		kappa = kappa + log(0.75)

		#survival probabilities
		surprob = read.table("~/git/abbg/R/old/KV/Input/surprobsmooth.txt")$V1 #60~89
		annprem <- surprob

		#unconditional survival prob
		unconsurprob = rep(1,Tret)
		for(it in 2:Tret){
		  unconsurprob[it] = unconsurprob[it-1]* surprob[it-1]
		}

    popsize <- rep(0, Ttot)
		for(it in 1:Twork){
    	popsize[it] = 1.0/(Twork + sum(unconsurprob))
		}

		for(it in 1:Tret){
		  popsize[Twork+it] = unconsurprob[it]/(Twork + sum(unconsurprob))
		}

		#####################################################
		#GRIDS
		#####################################################
		###################
		#Transitory shocks
		#lm <- optimize(FnGridTrans, c(0,4), p, tol=1e-4)$minimum
		ensd <- uniroot(FnGridTrans, c(1,4), p, extendInt="yes", tol=1e-2, maxiter=30)$root
		lval = FnGridTrans(ensd, p, FALSE) 
		edist <- lval$edist
		egrid <- lval$egrid 

		###################
		#Permanent component
		if (rho==1) {
    	Veta = Veta_rho1
		}else {  
    	Veta = (1-rho^2)*(Twork*Veta_rho1 - Vz0*(rho^(2*Twork)-1) ) /(1-rho^(2*Twork))
    }
		Vetavec  <- rep(Veta, Twork-1)
		
		varz   <- rep( Vz0, Twork	)
		for(it in 2:Twork){
		  varz[it] = (rho^2)*varz[it-1] + Vetavec[it-1]
		}
		
		#lval = comp.eta.prob(p,varz)
		#zdist      <- lval$zdist      
		#zgrid      <- lval$zgrid      
		#ztrans     <- lval$ztrans      
		#varzapprox <- lval$varzapprox 

		#save( zdist,zgrid,ztrans,varzapprox,file=paste('eta',ngpz,'dat',sep='.') )
		#load('~/git/abbg/R2/eta.dat')
    load( paste('~/git/abbg/R2/eta',ngpz,'dat',sep='.') )

		###################
		#Earnings
		#p$stax <- uniroot(FnTaxParamNet, c(0, 1), p, kappa, popsize, zgrid, egrid, zdist, edist,
    #  extendInt="yes", tol=1e-6, maxiter=200)$root
		#p$stax <-  0.003699068 
		p$stax <- 0.003704888
		lval <- FnTaxParamNet(p$stax, p, kappa, popsize, zgrid, egrid, zdist, edist, FALSE) 
		ygrid      <- lval$ygrid   
		ypregrid   <- lval$ypregrid
		avearnspre <- lval$avearnspre #avearnspost avlearnspre avlearnspost  avearnspre2 avearnspost2  
		  #avlearnspre2  avlearnspost2 varearnspre varearnspost  varlearnspre  varlearnspost 
	
		#Mean pre-tax earnings grid for pension
		p$pencap  <- pencapfrac* sum(avearnspre)/(Twork)

		###################
		#simulate to get distribution of average pre-tax incomes
		zsimI <- array( 0,dim=c(nsim, Twork) )
	  esimI <- zsimI
	  zsimI[,1] = sample( ngpz,nsim,replace=T,prob=zdist[1,] )
	  esimI[,1] = sample( ngpe,nsim,replace=T,prob=edist[1,] )

	  yavsim <- zsimI
	  ypresim <- array( 0,dim=c(nsim, Ttot) )

		ypresim[,1] = ypregrid[ cbind(1,zsimI[,1],esimI[,1]) ]
		yavsim[,1]  = ypresim[,1]
		yavsim[,1][ypresim[,1]>p$pencap] = p$pencap

		for (it in 2:Twork){
			for (i in 1:nsim){
		     zsimI[i,it]   = sample( ngpz,1,prob=ztrans[it-1,zsimI[i,it-1],] )
		  }	 
		  esimI[,it]   = sample( ngpe,nsim,replace=TRUE,prob=edist[it,] )
		  ypresim[,it] = ypregrid[ cbind(it,zsimI[,it],esimI[,it]) ]
		  yavsim[,it]  = ( (it-1)*yavsim[,it-1] + pmin(ypresim[,it],p$pencap) )/it
		}

		#equally spaced between 5th perc and median, and median and 95th perc
		mgrid    <- array( 0, dim=c(Twork,ngpm) ) #average pre-tax earnings grid
		vm <- seq(1/(2*ngpm), (2*ngpm-1)/(2*ngpm), l=ngpm) #median of each bin   
		for (it in 1:Twork) mgrid[it,] = quantile( yavsim[,it], vm, names=FALSE, na.rm = T )  

		###################
		#Pensions
		#sspar <- uniroot(FnSSParam, c(0.1,1.5), p, avearnspre, mgrid, extendInt="yes", tol=1e-3, maxiter=30)$root
		sspar <- 1.06032 
		lval <- FnSSParam(sspar, p, avearnspre, mgrid,FALSE)
		ppregrid <- lval$ppregrid
		pgrid    <- lval$pgrid  

		###################
		#Natural borrowing limits
		nbl	     <- rep( 0, Twork	)	# natural borrowing limits
		nblret	 <- array( 0, dim=c(Tret,ngpm) )	# natural borrowing limits

		#Last period
		for( im in 1:ngpm ) nblret[Tret,im] = -pgrid[Tret,im]/Rnet + 0.1

		#Retired
		for( it in (Tret-1):1 ){
		  for( im in 1:ngpm ){
		    nblret[it,im] = annprem[it]*nblret[it+1,im]/Rnet - pgrid[it,im]/Rnet + 0.1
		  }
		}

		#Last Working Period
		nbl[Twork] = nblret[1,1]/Rnet - ygrid[Twork,1,1]/Rnet + 0.1
		#Working
		for( it in (Twork-1):1 ) nbl[it] = nbl[it+1]/Rnet - ygrid[it,1,1]/Rnet +0.1

		###################
		#Assets: 0 to amax - exponentially spaced grid - pexpgrid is parametr
		agrid    <- array( 0, dim=c(Twork,ngpa)           ) #assets
 		agridret <- array( 0, dim=c(Tret,ngpa,ngpm)       ) #asset

 		for(it in 1:Twork){	
 			agrid[it,1] = max(borrowlim,nbl[it])
			for( ia in 2:ngpa ){
			  agrid[it,ia] = agrid[it,ia-1] + exp( pexpgrid*(ia-1) )
		  }
		  ltemp1 = ( amax-agrid[it,1] ) / ( agrid[it,ngpa]-agrid[it,1] )
		  agrid[it,] = agrid[it,1] + ltemp1*( agrid[it,]-agrid[it,1] )
		}  

    for(it in 1:Tret){  
      for( im in 1:ngpm ){
        agridret[it,1,im] = max( borrowlim,nblret[it,im] )
        for( ia in 2:ngpa ){
          agridret[it,ia,im] = agridret[it,ia-1,im] + exp( pexpgrid*(ia-1) )
        }
        ltemp1 = ( amax-agridret[it,1,im] ) / ( agridret[it,ngpa,im]-agridret[it,1,im] )
        agridret[it,,im] = agridret[it,1,im] + ltemp1*( agridret[it,,im]-agridret[it,1,im] )
      }
    } 

		###################
		#Cash in Hand
		xgrid    <- array( 0, dim=c(Twork,ngpa,ngpz,ngpe) ) #cash on hand
 		xgridret <- array( 0, dim=c(Tret,ngpa,ngpm)       )
 		tran <- xgrid
 		tranret  <- xgridret 

		for( ia in 1:ngpa ){
		  for( iz in 1:ngpz ){
		    for( ie in 1:ngpe ){
		      for( it in 1:(Twork-1) ){
		        tran[it,ia,iz,ie] = max(cfloor - (Rnet*agrid[it,ia] + ygrid[it,iz,ie] - agrid[it+1,1]), 0.0)
		        xgrid[it,ia,iz,ie] = Rnet*agrid[it,ia] + ygrid[it,iz,ie] + tran[it,ia,iz,ie]        
		      }
		      tran[Twork,ia,iz,ie] = max(cfloor - (Rnet*agrid[Twork,ia] + ygrid[Twork,iz,ie] - agridret[1,1,1]), 0.0)
		      xgrid[Twork,ia,iz,ie] = Rnet*agrid[Twork,ia] + ygrid[Twork,iz,ie] + tran[Twork,ia,iz,ie]        	      
		    }
		  }

      for( im in 1:ngpm ){
        for( it in 1:(Tret-1) ){
          tranret[it,ia,im] = max(cfloor - (Rnet*agridret[it,ia,im] + pgrid[it,im] - agridret[it+1,1,im] ),0.0 ) 
          xgridret[it,ia,im] = Rnet*agridret[it,ia,im] + pgrid[it,im] +tranret[it,ia,im]
        }
        tranret[Tret,ia,im] = max(cfloor - (Rnet*agridret[Tret,ia,im] + pgrid[Tret,im] ),0.0 ) 
        xgridret[Tret,ia,im] = Rnet*agridret[Tret,ia,im] + pgrid[Tret,im] +tranret[Tret,ia,im]
      }
		}

		if (Display==1) cat('Finished forming grids\n')

		#####################################################
		#Decisions: Retired
		#####################################################
		conret   <- xgridret
		mucret   <- xgridret
		assret   <- xgridret  #asset at the end of the period
		conret1  <- xgridret
		mucret1  <- xgridret 
		assret1  <- xgridret #asset at the beg
		          
		# Start with last period, eat everything  ##
		if (Display==1) cat( 'Solving for decision rules at age ', Ttot, '\n')
		mucret[Tret,,] = conret[Tret,,]^(-gam)           

		## Retired Agents: Pension value as state variable ##

		for( it in (Tret-1):1){
		  if (Display==1) cat( 'Solving for decision rules at age ', Twork+it, '\n')               

		  for( im in 1:ngpm ){          
		    #solve on tt+1 grid
		    mucret1[it,,im] = bet*(surprob[it]/annprem[it])*Rnet*mucret[it+1,,im]     
		    conret1[it,,im] = mucret1[it,,im]^(-1.0/gam)  

		    #asset level for desired consumption
		    assret1[it,,im] = (conret1[it,,im] +agridret[it+1,,im]*annprem[it]-pgrid[it,im] -tranret[it,,im])/Rnet #asset at the beg

		    #deal with borrowing limits
		    if ( min(agridret[it,,im]) >= assret1[it,1,im]) {     #my resource enough for desired consumption
		                                                          #BL does not bind anywhere
		      BLbind = 0
		    }else{                #find point in tt grid where BL starts to bind
		      BLbind = which.max( agridret[it,,im][ agridret[it,,im] < assret1[it,1,im] ] )
		      assret[it,1:BLbind,im]  = agridret[it+1,1,im] #spend as much as I can, end of period asset equals the min level of asset tom
		      conret[it,1:BLbind,im]  = xgridret[it,1:BLbind,im] - assret[it,1:BLbind,im]*annprem[it]
		      mucret[it,1:BLbind,im]  = conret[it,1:BLbind,im]^(-gam)
		    }
		    
		    #interpolate muc1 as fun of ass1 to get muc where BL does not bind
		    conret[it,(BLbind+1):ngpa,im] = approxExtrap( assret1[it,,im], conret1[it,,im], agridret[it,(BLbind+1):ngpa,im] )$y
		    mucret[it,(BLbind+1):ngpa,im] = conret[it,(BLbind+1):ngpa,im]^(-gam)       
		    assret[it,(BLbind+1):ngpa,im] = (xgridret[it,(BLbind+1):ngpa,im] - conret[it,(BLbind+1):ngpa,im])/annprem[it]

		    if ( any(!is.finite(conret[it,,im])) ) cat( 'nan encountered\n')
		  }
		}  

		#############################################           
		## Working Agents: Rules depend on shocks  ##
		#############################################
		con  <- array( 0, dim=c(Twork,ngpa,ngpm,ngpz,ngpe) )
		ass  <- con
		muc  <- con
		con1 <- con 
		ass1 <- con 
		muc1 <- con 

		model = list(	
			mucret   =mucret  ,  
			mgrid    =mgrid   , 
			ypregrid =ypregrid, 
			ygrid    =ygrid   ,    
			ztrans   =ztrans  ,  
			edist    =edist   , 
			agridret =agridret,    
			agrid    =agrid   , 
			tran     =tran    , 
			xgrid    =xgrid  
		)

		if(mode == 'mpi'){ 
			cat('[mode=mpi] USING MPI !!!!! \n')
			require(snow) 
			require(Rmpi) 
			cl <- makeCluster(type='MPI', spec=39)
			chainN = length(cl) 
			cat('Number of Chains: ',chainN,'\n')

			clusterEvalQ( cl,require(Hmisc) )
			clusterEvalQ( cl,source('~/git/abbg/R2/inc.modelsolver.nl.mpi.r') )
		} else if (mode == 'multicore'){
			cat('[mode=multicore] YEAH !!!!! \n')
			require(parallel)
			chainN = detectCores()
			cat('Number of Chains: ',chainN,'\n')
		}	

		for( it in Twork:1 ){
		  if (Display==1) cat( 'Solving for decision rules at age ', it, '\n')             

			if(mode == 'mpi'){ 
				vals <- parLapply(cl, 1:ngpz, comp.ngpz, p, model, muc, it)
		  }else if (mode == 'multicore'){
        vals <- mclapply(1:ngpz, comp.ngpz, p, model, muc, it, mc.cores = chainN )
      }else{
        vals <- lapply(1:ngpz, comp.ngpz, p, model, muc, it)
      } 

			for(iz in 1:ngpz){
				 con[it, , ,iz, ] <- vals[[iz]]$lcon  
				 ass[it, , ,iz, ] <- vals[[iz]]$lass  
				 muc[it, , ,iz, ] <- vals[[iz]]$lmuc  
				con1[it, , ,iz, ] <- vals[[iz]]$lcon1
				ass1[it, , ,iz, ] <- vals[[iz]]$lass1
				muc1[it, , ,iz, ] <- vals[[iz]]$lmuc1
			}
		}  #t

		if (mode == 'mpi') stopCluster(cl)

		models= list(
			popsize = popsize,
			annprem = annprem, 

			ztrans  = ztrans ,
			zgrid   = zgrid  ,
			zdist   = zdist  ,

			egrid   = egrid  ,
			edist   = edist  ,  

			ypregrid= ypregrid,
			ygrid   = ygrid   ,

			mgrid   = mgrid  ,
			
			agrid   = agrid  , 
			agridret= agridret, 
			
			pgrid   = pgrid  ,
			ppregrid= ppregrid,

			con     = con    ,
			conret  = conret ,
			p       = p
		)   
	}) 
  return(res)
}

#####################################################
#Simulations until age 35, and create sample
#####################################################
comp.samples <- function(models, tau0){
	res <- with( c(models, models$p), { 
		
		zsimI  <- array( 0,dim=c(nsim, Twork) )
	  esimI  <- zsimI
	  yavsim <- zsimI
		msimI  <- zsimI
		esim   <- zsimI
		zsim   <- zsimI
		msim   <- zsimI

	  ypresim <- array( 0,dim=c(nsim, Ttot) )
		ysim    <- ypresim
		csim    <- ypresim
		xsim    <- ypresim
		trsim   <- ypresim

		asim    <- array( 0,dim=c(nsim, Ttot+1) )

		#initial earnings
	  zsimI[,1] = sample( ngpz,nsim,replace=T,prob=zdist[1,] )
	  zsim[,1]  = zgrid[1,zsimI[,1]]

	  esimI[,1] = sample( ngpe,nsim,replace=T,prob=edist[1,] )
		esim[,1]  = egrid[1,esimI[,1]]

		ysim[,1]    = ygrid[ cbind(1,zsimI[,1],esimI[,1]) ]
		ypresim[,1] = ypregrid[ cbind(1,zsimI[,1],esimI[,1]) ]
		yavsim[,1]  = ypresim[,1]
		yavsim[,1][ypresim[,1] > pencap] = pencap

		for( i in 1:nsim ){
		  im2 = FindLinProb1(yavsim[i,1],mgrid[1,])
		  msimI[i,1] = sample( c(im2[1],im2[1]+1), 1, prob=im2[2:3] )
		}
		msim[,1] = mgrid[1,msimI[,1]]

		#do t=1 separately
		it = 1
		for( i in 1:nsim ){
		  trsim[i,it] = max( cfloor - (Rnet*asim[i,it]+ ysim[i,it] - agrid[it+1,1]) ,0.0 )
		  csim[i,it]  = approxExtrap( agrid[it, ], con[ it,, msimI[i,it], zsimI[i,it], esimI[i,it] ], asim[i,it]+trsim[i,it] )$y
		  if (trsim[i,it]>0.0) csim[i,it] = min(cfloor, csim[i,it])
		  xsim[i,it]     = Rnet*asim[i,it]+ ysim[i,it] + trsim[i,it]
		  asim[i,it+1]   = xsim[i,it] - csim[i,it]
		  if ( asim[i,it+1]<agrid[it+1,1] ){
		    asim[i,it+1] = agrid[it+1,1]
		    csim[i,it] = xsim[i,it] - asim[i,it+1]
		  }        
		  if(!is.finite(csim[i,it]))  cat('bad consumption\n')    
		}

		#working life
		for( it in 2:irb ){
			for (i in 1:nsim) zsimI[i,it] = sample( ngpz,1,prob=ztrans[it-1,zsimI[i,it-1],] )
		  zsim[,it] = zgrid[it,zsimI[,it]]

		  esimI[,it]   = sample( ngpe,nsim,replace=TRUE,prob=edist[it,] )
		  esim[,it] = egrid[it,esimI[,it]]

		  ysim[,it] = ygrid[ cbind(it,zsimI[,it],esimI[,it]) ]
		  ypresim[,it] = ypregrid[ cbind(it,zsimI[,it],esimI[,it]) ]
		  yavsim[,it]  = ( (it-1)*yavsim[,it-1] + pmin(ypresim[,it], pencap) )/it

		  for( i in 1:nsim ){
		    im2 = FindLinProb1(yavsim[i,it], mgrid[it,])
		    msimI[i,it] = sample( c(im2[1],im2[1]+1), 1, prob=im2[2:3] )
		    msim[i,it] = mgrid[it, msimI[i,it]]
		      
		    if (it<Twork) {
		      trsim[i,it]    = max( cfloor - (Rnet*asim[i,it]+ ysim[i,it] - agrid[it+1,1]) ,0.0)
		    }else if(it==Twork) {
		      trsim[i,it]    = max( cfloor - (Rnet*asim[i,it]+ ysim[i,it] - agridret[1,1,msimI[i,Twork]] ) ,0.0)
		    }
		    csim[i,it]  = approxExtrap( agrid[it, ], con[ it,, msimI[i,it], zsimI[i,it], esimI[i,it] ], asim[i,it] )$y
		    xsim[i,it]     = Rnet*asim[i,it]+ ysim[i,it] + trsim[i,it]
		    asim[i,it+1]   = xsim[i,it] - csim[i,it]
		    
		    if (it<Twork) {    
		      if( asim[i,it+1] < agrid[it+1,1] ){
		        asim[i,it+1] = agrid[it+1,1]
		        csim[i,it] = xsim[i,it] - asim[i,it+1]
		      }        
		    }else if(it==Twork) {
		      if( asim[i,it+1] < agridret[ 1,1,msimI[i,Twork] ] ){
		        asim[i,it+1] = agridret[ 1,1,msimI[i,Twork] ]
		        csim[i,it] = xsim[i,it] - asim[i,it+1]
		      }                
		    }
		    
		    if(!is.finite(csim[i,it])) cat('bad consumption\n')
		  }
		}

		# Keep only tau0-percentile of initial eta
		eta_1 = quantile(zsim[,irb], tau0)
		hhnum <- which(zsim[,irb] == eta_1)
		if( length(hhnum) < 400 ) cat('need more in percentile ',tau0,'\n')
		
		indnum <- sample(hhnum, 400) #random select 1000 individuals
		hhnum <- rep(indnum, each = 125) #enlarge to original sample

		#replace values before 35
		zsimI[,1:irb] <- zsimI[hhnum,1:irb]
		zsim[,1:irb]  <- zsim[hhnum,1:irb]

		esimI[,1:irb] <- esimI[hhnum,1:irb]
		esim[,1:irb]  <- esim[hhnum,1:irb]
				
		msimI[,1:irb] <- msimI[hhnum,1:irb]
		msim[,1:irb]  <- msim[hhnum,1:irb]

		ypresim[,1:irb] <- ypresim[hhnum,1:irb]
		yavsim[,1:irb]  <- yavsim[hhnum,1:irb]
		ysim[,1:irb]    <- ysim[hhnum,1:irb]

		csim[,1:irb]     <- csim[hhnum,1:irb]
		xsim[,1:irb]     <- xsim[hhnum,1:irb]
		asim[,1:(irb+1)] <- asim[hhnum,1:(irb+1)]
		trsim[,1:irb]    <- trsim[hhnum,1:irb]

		samples = list(
		  zsimI   = zsimI  ,
		  zsim    = zsim   , 

		  esimI   = esimI  , 
			esim    = esim   , 
		  
		  msimI   = msimI  , 
		  msim    = msim   , 
		  
		  yavsim  = yavsim , 
			ysim    = ysim   , 
			ypresim = ypresim,
			
			asim    = asim   , 
			xsim    = xsim   , 
			csim    = csim   ,
			trsim   = trsim
		)   
	}) 
  return(res)
}

#####################################################
#compute moments
#####################################################
comp.moments <- function(models, samples, tau0, tau1) {
	res <- with( c(models, models$p, samples), { 

		#find shock tau1=0.5 at age 36
		load('~/git/abbg/R2/sig_v.dat')
    shockperm <-  qnorm( tau1, sd=sig_v[irb] )
    iab  = (tau0 == 0.9) & (tau1 == 0.1)
    ibe  = (tau0 == 0.1) & (tau1 == 0.9)
    iyes = iab | ibe    

    zshock <- (1-delta*iyes) * unique(zsim[,irb]) + shockperm
		iz2 <- FindLinProb1( zshock, zgrid[irb+1,] )

		#working life after shock
		for ( it in (irb+1):Twork ){
			if( it==(irb+1) ) {
				zsimI[,it] = sample( c(iz2[1],iz2[1]+1), nsim, replace=TRUE, prob=iz2[2:3] )
			}else{
				for (i in 1:nsim){
			     zsimI[i,it]   = sample( ngpz,1,prob=ztrans[it-1,zsimI[i,it-1],] )
			  }
			}  	 
			zsim[,it]    = zgrid[it,zsimI[,it]]

		  esimI[,it]   = sample( ngpe,nsim,replace=TRUE,prob=edist[it,] )
		  esim[,it]    = egrid[it,esimI[,it]]

		  ysim[,it]    = ygrid[ cbind(it,zsimI[,it],esimI[,it]) ]
			ypresim[,it] = ypregrid[ cbind(it,zsimI[,it],esimI[,it]) ]
			yavsim[,it]  = ( (it-1)*yavsim[,it-1] + pmin(ypresim[,it], pencap) )/it

		  for( i in 1:nsim ){
		    im2 = FindLinProb1(yavsim[i,it], mgrid[it,])
		    msimI[i,it] = sample( c(im2[1],im2[1]+1), 1, prob=im2[2:3] )
		    msim[i,it] = mgrid[it, msimI[i,it]]
		      
		    if (it<Twork) {
		      trsim[i,it]    = max( cfloor - (Rnet*asim[i,it]+ ysim[i,it] - agrid[it+1,1]) ,0.0)
		    }else if(it==Twork) {
		      trsim[i,it]    = max( cfloor - (Rnet*asim[i,it]+ ysim[i,it] - agridret[1,1,msimI[i,Twork]] ) ,0.0)
		    }
		    csim[i,it]  = approxExtrap( agrid[it, ], con[ it,, msimI[i,it], zsimI[i,it], esimI[i,it] ], asim[i,it] )$y
		    xsim[i,it]     = Rnet*asim[i,it]+ ysim[i,it] + trsim[i,it]
		    asim[i,it+1]   = xsim[i,it] - csim[i,it]
		    
		    if (it<Twork) {    
		      if( asim[i,it+1] < agrid[it+1,1] ){
		        asim[i,it+1] = agrid[it+1,1]
		        csim[i,it] = xsim[i,it] - asim[i,it+1]
		      }        
		    }else if(it==Twork) {
		      if( asim[i,it+1] < agridret[ 1,1,msimI[i,Twork] ] ){
		        asim[i,it+1] = agridret[ 1,1,msimI[i,Twork] ]
		        csim[i,it] = xsim[i,it] - asim[i,it+1]
		      }                
		    }
		    
		    if(!is.finite(csim[i,it])) cat('bad consumption\n')
		  }
		}

		#retirement
		for( it in 1:Tret ){
		  for( i in 1:nsim ){
		    ysim[i,Twork+it] = pgrid[ it,msimI[i,Twork] ]
		    ypresim[i,Twork+it] = ppregrid[ it,msimI[i,Twork] ]
		    if (it<Tret) trsim[i,Twork +it]    = max( cfloor - ( Rnet*asim[i,it]+ ysim[i,it] - agridret[ it+1,1,msimI[i,Twork] ] ) ,0.0 )        
		    csim[i,Twork+it] = approxExtrap( agridret[ it,, msimI[i,Twork] ], conret[ it, ,msimI[i,Twork] ], asim[i,Twork+it] )$y
		    xsim[i,Twork+it]   = Rnet*asim[i,Twork+it]+ ysim[i,Twork+it]  + trsim[i,Twork +it]
		    
		    if (it<Tret) {
		      asim[i,Twork+it+1] = (xsim[i,Twork+it]  - csim[i,Twork+it])/annprem[it]
		      if ( asim[i,Twork+it+1] < agridret[ it+1,1,msimI[i,Twork] ] ){
		        asim[i,Twork+it+1] = agridret[ it+1,1,msimI[i,Twork] ]
		        csim[i,Twork+it] = xsim[i,Twork+it] - asim[i,Twork+it+1]*annprem[it]
		      }        
		    }else if(it==Tret) {
		      asim[i,Twork+it+1] = xsim[i,Twork+it]  - csim[i,Twork+it]
		      if (asim[i,Twork+it+1] != 0){
		        asim[i,Twork+it+1] = 0
		        csim[i,Twork+it] = xsim[i,Twork+it]
		      }             
		    }
		  }
		}

		moments= list(
		  zsim    = zsim   , 
			esim    = esim   , 

			ysim    = ysim   , 
			ypresim = ypresim,
			
			asim    = asim   , 
			csim    = csim
		)   
	}) 

  return(res)
}

