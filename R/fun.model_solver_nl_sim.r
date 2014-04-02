source('inc.modelsolver_nl.r')

comp.eta <- function(p){
	res <- with(p,{  
		# Number of individuals
		N=1000000
		#N=1000
		aa_ref=30
		V_draw <- runif(N)  

		MatAGE1 <- rep(0, K3+1)
		for (kk3 in 0:K3){    
		  MatAGE1[kk3+1]=hermite( (aa_ref-meanAGE)/stdAGE, kk3 )
		}
		Mateta_true = array(0, dim = c(N,nage))

		Mateta_true[,1] = (MatAGE1 %*% Resqtrue_e0[,1]) * (V_draw <= Vectau[1])
		for (jtau in 2:Ntau){
	    Mateta_true[,1] = Mateta_true[,1] + ( (MatAGE1 %*% 
	    	(Resqtrue_e0[,jtau]-Resqtrue_e0[,jtau-1]))/(Vectau[jtau]-Vectau[jtau-1]) *
	      (V_draw-Vectau[jtau-1]) + MatAGE1 %*% Resqtrue_e0[,jtau-1] ) * 
	      (V_draw>Vectau[jtau-1]) * (V_draw<=Vectau[jtau])
		}
		Mateta_true[,1]=Mateta_true[,1]+(MatAGE1 %*% Resqtrue_e0[,Ntau])*(V_draw>Vectau[Ntau])

		Mateta_true[,1]=Mateta_true[,1]+( (1/(b1true_e0)*log(V_draw/Vectau[1]))*(V_draw<=Vectau[1]) - 
		  (1/bLtrue_e0*log((1-V_draw)/(1-Vectau[Ntau])))*(V_draw>Vectau[Ntau]))

		for (jj in 1:nage){
		    
		  aa=aa_ref+(jj-1)*2
		   
	    # Eta		    
	    if (jj <= nage-1){ #this period is jj+1
	        
        Mat = array( 0, dim=c(N, (K1+1)*(K2+1)) )
        for (kk1 in 0:K1){
          for (kk2 in 0:K2) {
            Mat[,kk1*(K2+1)+kk2+1] = hermite( (Mateta_true[,jj]-meanY)/stdY, kk1 ) * 
            	hermite( ((aa+2)-meanAGE)/stdAGE, kk2 )
          }
        }
	        
        V_draw <- runif(N)  
        
        #First quantile
        
        Mateta_true[,jj+1] = ( Mat %*% Resqtrue[,1] ) * ( V_draw <= Vectau[1] )

        for (jtau in 2:Ntau){
          Mateta_true[,jj+1] = Mateta_true[,jj+1] + 
            ( (Mat %*% (Resqtrue[,jtau]-Resqtrue[,jtau-1])) / 
            (Vectau[jtau]-Vectau[jtau-1]) *
            (V_draw - Vectau[jtau-1]) + Mat %*% Resqtrue[,jtau-1] ) * 
            (V_draw > Vectau[jtau-1]) * (V_draw <= Vectau[jtau])
        }
        Mateta_true[,jj+1]=Mateta_true[,jj+1] + (Mat %*% Resqtrue[,Ntau]) * (V_draw>Vectau[Ntau])

        Mateta_true[,jj+1]=Mateta_true[,jj+1] + ( (1 / b1true * log(V_draw/Vectau[1])) * (V_draw <= Vectau[1]) - 
            (1 / bLtrue * log((1-V_draw)/(1-Vectau[Ntau]))) * (V_draw > Vectau[Ntau]) )
	    }
		}

		# Quantiles of eta and epsilon, by age
		xeta <- array( 0, dim=c(nage, nbin) )  #23 bin
		veta <- seq(1/(2*nbin), (2*nbin-1)/(2*nbin), l=(2*nbin-1)) #45 nodes
		oddnode <- seq(1,(2*nbin-1),2)  #income node is the even node
		age = seq(30, (30+(nage-1)*2), 2)

		for (i in 1:nage){
		  xeta[i,] <- quantile( Mateta_true[,i], veta, names=F, na.rm = T )[oddnode]
		}    

		long <- data.table( pid = 1:N, age=rep(age, each=N), income = c(Mateta_true))
		setkey(long, pid, age)
		long[, q:=as.numeric(cut_number(income, n = nbin)), age]
		long[,income:= NULL]
		wide <- reshape(long, idvar='pid', timevar='age', direction='wide')

		etaprob <- array( 0,dim=c(nage, nbin, nbin) )
		etaprob[1,,] <- 1/nbin

		for (t in 1:(nage-1)){
		  x1 <- paste('q',age[t],sep='.')
		  x2 <- paste('q',age[t]+2,sep='.')
		  etaprob[t+1,,] <- trans.matrix( wide[[x1]], wide[[x2]] )
		}    

		etacontot <- etaprob
		for (i in 1:nage) {
		  etacontot[i,,] <- econdCDF( nbin, etaprob[i,,] )
		}  

		etauntot <- etaprob[1,1,]
		for(i in 2:nbin){
		  etauntot[i] <- etauntot[i-1]  + etauntot[i]  
		}

		list(
		xeta      = xeta,
		etaprob   = etaprob,
    etacontot = etacontot,
    etauntot  = etauntot)
	
	}) 
  return(res)
}


comp.solveModel <- function(p) {
	res <- with(p,{  
		
		#Setting up shock values (discrete approximation to log normal)   
		PermVecProb = rep(1/nP,nP)
		ThetaVecProb = rep(1/nT,nT)
		if (pUnemp>0) { # If assume unemployment
   	 ThetaVecProb = c( pUnemp,ThetaVecProb*(1-pUnemp) )
 		}

		PermVec  = DiscreteApproxToMeanOneLogNormal(sigP,nP)
		ThetaVec = DiscreteApproxToMeanOneLogNormal(sigT,nT)
		if (pUnemp>0){ # If assume unemployment
    	ThetaVec = c( 0,ThetaVec/(1-pUnemp) )
 		}

		ThetaMat = t( replicate(40, ThetaVec) )
		ThetaMat = rbind( ThetaMat, matrix(1,nrow=25,ncol=length(ThetaVec)) )

		PermMat = t( replicate(40, PermVec) )
		PermMat = rbind( PermMat, matrix(1,nrow=25,ncol=length(PermVec)) )

    AlphaVec=exp(exp(exp(seq(log(log(log(aMin+1)+1)+1),log(log(log(aMax+1)+1)+1),l=n))-1)-1)-1  # set up triple exponential grid
		AlphaVec = c(AlphaVec,aHuge)

		## Profile of income level 
		Incomec = c(1.0304073e+001,1.0343448e+001,1.0382479e+001,1.0421080e+001,1.0459171e+001,1.0496674e+001,1.0533513e+001,1.0569617e+001,1.0604917e+001,1.0639347e+001,1.0672845e+001,1.0705351e+001,1.0736808e+001,1.0767163e+001,1.0796365e+001,1.0824368e+001,1.0851126e+001,1.0876600e+001,1.0900750e+001,1.0923542e+001,1.0944944e+001,1.0964928e+001,1.0983466e+001,1.1000538e+001,1.1016122e+001,1.1030203e+001,1.1042768e+001,1.1053805e+001,1.1063308e+001,1.1071273e+001,1.1077697e+001,1.1082585e+001,1.1085939e+001,1.1087769e+001,1.1088086e+001,1.1086903e+001,1.1084239e+001,1.1080114e+001,1.1074551e+001,1.1067577e+001,1.1059222e+001,1.0369837e+001,1.0369227e+001,1.0368617e+001,1.0368007e+001,1.0367397e+001,1.0366786e+001,1.0366176e+001,1.0365566e+001,1.0364956e+001,1.0364345e+001,1.0363735e+001,1.0363125e+001,1.0362515e+001,1.0361904e+001,1.0361294e+001,1.0360684e+001,1.0360074e+001,1.0359463e+001,1.0358853e+001,1.0358243e+001,1.0357633e+001,1.0357023e+001,1.0356412e+001,1.0355802e+001,1.0355192e+001)
		IncomeLevelc = exp(Incomec)

		# Income growth factors
		GList =  IncomeLevelc[-1]/IncomeLevelc[-length(IncomeLevelc)]
	

		# Construct matrix of interpolation data
		C = matrix(0:(n+1))
		M = matrix(0:(n+1))

		for (l in 1:PeriodsToSolve){			
			# calculate ct from each grid point in AlphaVec
			Vap    = GothVP(p, AlphaVec, GList[65-l+1], ThetaMat[65-l+1,], PermMat[65-l+1,], ThetaVecProb, PermVecProb, C, M )   # Gothic Va prime
			ChiVec = Vap^(-1/p$rho) # inverse Euler equation
 			MuVec  = AlphaVec+ChiVec
  		M      = cbind(M, c(0,MuVec))                  # Matrix of interpolation data
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

		# Construct income shock draw lists
		# Tran shock draw list
		ThetaDraws = DiscreteApproxToMeanOneLogNormal(sigT,NumOfPeople)
		if (pUnemp>0){ # If assume unemployment
			ThetaDraws = DiscreteApproxToMeanOneLogNormal(sigT,NumOfPeople-pUnemp*NumOfPeople)
			ThetaDraws = ThetaDraws/(1 - pUnemp)
			ThetaDraws = append(rep(0,pUnemp*NumOfPeople),ThetaDraws)
		} 

		# Perm shock draw list
		PermShockDraws = DiscreteApproxToMeanOneLogNormal(sigP,NumOfPeople)

		# Construct income shock lists
		for (i in 1:40){
		  ThetaList[i,] = sample(ThetaDraws)     # List of Theta (tran shock) 
		  PermList[i,]  = sample(PermShockDraws) # List of perm shock 
		}
		Perm[1,] = PermList[1,]
		PermList[41:66,] = matrix(1,nrow=26,ncol=length(PermShockDraws)) 
		ThetaList[41:66,]= matrix(1,nrow=26,ncol=length(ThetaDraws)) 

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

		len = ncol(M)
		ctList[1,] =  approx( M[,len], C[,len], mtList[1,] )$y
		ctList[1, mtList[1,] == 0] =0                     

		# Continue simulstion
		for (t in 2:NumOfPeriodsToSimulate){
	    for (j in 1:NumOfPeople){           
        stList[t,j] = ( R/GList[t-1]/PermList[t,j] )*( mtList[t-1,j]-ctList[t-1,j] )
        mtList[t,j] = stList[t,j] + ThetaList[t,j]
        ctList[t,j] = approx( M[,len-t+1], C[,len-t+1], mtList[t,j] )$y
        if (mtList[t,j] == 0) ctList[t,j] =0                       # ctList      : list of normalized consumption 
				Perm[t,j]   = Perm[t-1,j] * GList[t-1] * PermList[t,j] # List of perm income
			}
		}

		stMeanList   = rowMeans(stList)         # stLevelMedianList: list of median of savings level 
		stMedianList = apply(stList, 1, median)
		ctMeanList   = rowMeans(ctList) 
		ctMedianList = apply(ctList, 1, median)
		Ct           = ctList * Perm
		CtMean       = rowMeans(Ct)
		
		model= list(
		CtMean        = CtMean,
    stMeanList    = stMeanList,
    stMedianList  = stMedianList,
    ctMeanList    = ctMeanList,
    ctMedianList  = ctMedianList)
	}) 

  return(res)
}