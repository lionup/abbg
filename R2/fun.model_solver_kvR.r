source('inc.modelsolver_kvR.r')

comp.solveModel <- function(p) {
	res <- with(p,{  

		#####################################################
		#DECLARE
		#####################################################
		#FOR GRIDS
		agrid    <- array( 0, dim=c(Twork,ngpa)           ) #assets
   
		mgrid    <- array( 0, dim=c(Twork,ngpm)           ) #average pre-tax earnings grid   
		xgrid    <- array( 0, dim=c(Twork,ngpa,ngpz,ngpe) ) #cash on hand

		pgrid    <- array( 0, dim=c(Tret,ngpp)            ) #pension grid
		ppregrid <- pgrid
 		agridret <- array( 0, dim=c(Tret,ngpa,ngpp)       ) #assets
 		xgridret <- agridret

 		nbl	     <- rep( 0, Twork	)	# natural borrowing limits
		nblret	 <- array( 0, dim=c(Tret,ngpp)            )	# natural borrowing limits

		#GLOBALS TO STORE DECISION RULES 
		con <- array( 0, dim=c(Twork,ngpa,ngpm,ngpz,ngpe) )
		ass  <- con
		ass1 <- con 
		con1 <- con 
		muc  <- con
		muc1 <- con 
		emuc <- con 
		val  <- con
		tran <- xgrid
		conret   <- xgridret
		assret   <- xgridret
		conret1  <- xgridret 
		assret1  <- xgridret 
		mucret   <- xgridret
		mucret1  <- xgridret 
		Emucret1 <- xgridret  
		valret   <- xgridret
		tranret  <- xgridret 

		#GLOBALS FOR PENSION SYSTEM
    pind <- array( 0, dim=c(ngpm,ngpz,ngpe) )

    #var of permanet shock
    if (rho==1) {
    	Veta = Veta_rho1
		}else {  
    	Veta = (1-rho^2)*(Twork*Veta_rho1 - Vz0*(rho^(2*Twork)-1) ) /(1-rho^(2*Twork))
    }
		Vetavec  <- rep(Veta, Twork-1)

		#####################################################
		#DATA
		#####################################################
		#log experience profile
		kappa = read.table("../R/old/KV/Input/kappasmooth.txt")$V1 #25~59
		kappa = kappa + log(0.75)

		#survival probabilities
		surprob = read.table("../R/old/KV/Input/surprobsmooth.txt")$V1 #60~89
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
		#Transitory shocks
		lm <- optimize(FnGridTrans, c(0,4), p, tol=1e-4)$minimum
		elist = FnGridTrans(lm, p, FALSE) #edist, egrid 

		#Permanent component
		varz   <- nbl
		varz[1] = Vz0
		for(it in 2:Twork){
		  varz[it] = (rho^2)*varz[it-1] + Vetavec[it-1]
		}
		
		#start_time = proc.time()[3]  
		lm = optimize(FnGridPerm, c(1,4), p, varz, tol=1e-2)$minimum
		#cat(paste('\ntotal seconds to compute moments' , proc.time()[3] -  start_time ))
		zlist = FnGridPerm(lm, p, varz, FALSE) #zdist zgrid ztrans varzapprox 

		#Earnings
		stax <- uniroot(FnGrossInc, c(0, 1), p, kappa, popsize, zlist$zgrid, elist$egrid, zlist$zdist, elist$edist,
      extendInt="yes", tol=1e-6, maxiter=200)$root

		lval <- FnGrossInc(stax, p, kappa, popsize, zlist$zgrid, elist$egrid, zlist$zdist, elist$edist, FALSE) #ygrid, ypregrid 
		  #avearnspre, avearnspost avlearnspre avlearnspost  avearnspre2 avearnspost2  
		  #avlearnspre2  avlearnspost2 varearnspre varearnspost  varlearnspre  varlearnspost 
	
		#Mean pre-tax earnings grid for pension
		pencap = pencapfrac* sum(avearnspre)/(Twork)

		#simulate to get distribution of average pre-tax incomes

		#equally spaced between 5th perc and median, and median and 95th perc

		#Pensions


		
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

		#GLOBALS TO STORE SIMULATION RESULTS
		ysim    <- array( 0, dim=c(nsim, Ttot) )
		ypresim <- ysim   
		csim    <- ysim
		xsim    <- ysim
		trsim   <- ysim 
		esim    <- array( 0, dim=c(nsim, Twork) )
		zsim    <- esim
		tsim    <- esim
		yavsim  <- esim  
		msim    <- esim
		asim    <- array( 0, dim=c(nsim, Ttot+1) )
		esimI   <- esim
		zsimI   <- esim
		msimI   <- esim
		psimI   <- rep(0, nsim)

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