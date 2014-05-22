source('inc.modelsolver_nl_st.r')

comp.eta.sim <- function(p){
  res <- with(p,{  

    aa_ref=age_min
    #V_draw <- runif(N)  
    Vgrid<- (1:N) / (1+N)
    V_draw <- array(0, dim = c(nage,N))

    for (i in 1:nage){
      V_draw[i,] = sample(Vgrid)
    }

    MatAGE1 <- rep(0, K3+1)
    for (kk3 in 0:K3){    
      MatAGE1[kk3+1]=hermite( (aa_ref-meanAGE)/stdAGE, kk3 )
    }
    Mateta_true = array(0, dim = c(nage,N))

    Mateta_true[1,] = (MatAGE1 %*% Resqtrue_e0[,1]) * (V_draw[1,] <= Vectau[1])
    for (jtau in 2:Ntau){
      Mateta_true[1,] = Mateta_true[1,] + ( (MatAGE1 %*% 
        (Resqtrue_e0[,jtau]-Resqtrue_e0[,jtau-1]))/(Vectau[jtau]-Vectau[jtau-1]) *
        (V_draw[1,]-Vectau[jtau-1]) + MatAGE1 %*% Resqtrue_e0[,jtau-1] ) * 
        (V_draw[1,]>Vectau[jtau-1]) * (V_draw[1,]<=Vectau[jtau])
    }
    Mateta_true[1,]=Mateta_true[1,]+(MatAGE1 %*% Resqtrue_e0[,Ntau])*(V_draw[1,] > Vectau[Ntau])

    Mateta_true[1,]=Mateta_true[1,]+( (1/(b1true_e0)*log(V_draw[1,]/Vectau[1]))*(V_draw[1,]<=Vectau[1]) - 
      (1/bLtrue_e0*log((1-V_draw[1,])/(1-Vectau[Ntau])))*(V_draw[1,]>Vectau[Ntau]))

    for (jj in 1:nage){  #periods before retirement
        
      aa=aa_ref+(jj-1)*2 #age last period
       
      # Eta       
      if (jj <= nage-1){ #already have eta at period 1, need eta for the rest nage-1 period
                         #this period is jj+1
          
        Mat = array( 0, dim=c(N, (K1+1)*(K2+1)) )
        for (kk1 in 0:K1){
          for (kk2 in 0:K2) {
            Mat[,kk1*(K2+1)+kk2+1] = hermite( (Mateta_true[jj,]-meanY)/stdY, kk1 ) * 
              hermite( ((aa+2)-meanAGE)/stdAGE, kk2 )
          }
        }
          
        #V_draw <- runif(N) 
        #First quantile
        Mateta_true[jj+1,] = ( Mat %*% Resqtrue[,1] ) * ( V_draw[jj+1,] <= Vectau[1] )

        for (jtau in 2:Ntau){
          Mateta_true[jj+1,] = Mateta_true[jj+1,] + 
            ( (Mat %*% (Resqtrue[,jtau]-Resqtrue[,jtau-1])) / 
            (Vectau[jtau]-Vectau[jtau-1]) *
            (V_draw[jj+1,] - Vectau[jtau-1]) + Mat %*% Resqtrue[,jtau-1] ) * 
            (V_draw[jj+1,] > Vectau[jtau-1]) * (V_draw[jj+1,] <= Vectau[jtau])
        }
        Mateta_true[jj+1,]=Mateta_true[jj+1,] + (Mat %*% Resqtrue[,Ntau]) * (V_draw[jj+1,]>Vectau[Ntau])

        Mateta_true[jj+1,]=Mateta_true[jj+1,] + ( (1 / b1true * log(V_draw[jj+1,]/Vectau[1])) * (V_draw[jj+1,] <= Vectau[1]) - 
            (1 / bLtrue * log((1-V_draw[jj+1,])/(1-Vectau[Ntau]))) * (V_draw[jj+1,] > Vectau[Ntau]) )
      }
    }

    save_Mateta_name <- paste('Mateta',age_min,'.dat',sep='')
    save(Mateta_true, file=save_Mateta_name)
    
    Mateta_true
  }) 
  return(res)
}

comp.eta.prob <- function(p){
  res <- with(p,{ 

    # get the simulations of workers
    Mateta_true <- comp.eta.sim(p)
    save_Mateta_name <- paste('Mateta',age_min,'.dat',sep='')
    load(save_Mateta_name)

    # Quantiles of eta and epsilon, by age
    xeta <- array( 0, dim=c(nage, nbin) )  #bins
    veta <- seq(1/(2*nbin), (2*nbin-1)/(2*nbin), l=(2*nbin-1)) #bins and its median
    oddnode <- seq(1,(2*nbin-1),2)  #income node is the median in each interval, odd node
    age = seq(age_min, age_max, 2)

    for (i in 1:nage){ #periods before retirement
      xeta[i,] <- quantile( Mateta_true[i,], veta, names=F, na.rm = T )[oddnode]
    }  

    long <- data.table( pid = 1:N, age=rep(age, each=N), income = c(t(Mateta_true)) )
    setkey(long, pid, age)
    long[, q:=as.numeric(cut_number(income, n = nbin)), age] #give a bin # for each person each age
    long$income <- NULL
    wide <- reshape(long, idvar='pid', timevar='age', direction='wide') #each person has all age in a row

    etaprob <- array( 0,dim=c(nage, nbin, nbin) ) #last period to this period
    etaprob[1,,] <- 1/nbin  #first period, same prob

    for (t in 2:nage){ #today
      x1 <- paste('q',age[t-1],sep='.')
      x2 <- paste('q',age[t],sep='.')
      etaprob[t,,] <- trans.matrix( wide[[x1]], wide[[x2]] ) #today
    }     

    etacontot <- etaprob
    for ( i in 1:nage ) {
      etacontot[i,,] <- econdCDF( nbin, etaprob[i,,] )
    }  

    etauntot <- etaprob[1,1,]
    for(i in 2:nbin){
      etauntot[i] <- etauntot[i-1]  + etauntot[i]  
    }
    
    ieta = exp(xeta)

    list(
    xeta      = xeta,
    ieta      = ieta,
    etaprob   = etaprob,
    etacontot = etacontot,
    etauntot  = etauntot)
  
  }) 
  return(res)
}

comp.eps <- function(p, ntra){ 
  xeps <- array( 0, dim=c(p$nage, ntra) ) #income grid
  VecTau  <- (1:p$Ntau)/(1+p$Ntau) #get the quantile of interpolation node
  min <- 1/(1+p$neps)
  max <- p$neps/(1+p$neps)
  VecTaue <- seq(min, max, l=ntra)
  #VecTaue <- (1:ntra) / (1+ntra) #get the quantile of interpolation node

  # generate the age vector 
  xt  <- seq( p$age_min, (p$age_min + (p$nage-1) * 2 ), by=2 )
  xtL <- ( xt - data$meanAGE ) / data$stdAGE #standardize AGE

  # Create MatAGE which is the hermite of AGE
  MatAGE <- array( 0,dim=c(p$K4+1, p$nage) )
  for (kk4 in 0:p$K4){    
    MatAGE[kk4+1,]=hermite(xtL,kk4)
  }

  for (l in 1:p$nage){
    xeps[l,] = ( t(MatAGE[,l]) %*% p$Resqtrue_eps[,1] )*( VecTaue <= VecTau[1] )
    for (jtau in 2:p$Ntau){
        xeps[l,] = xeps[l,] + ( (t(MatAGE[,l]) %*% (p$Resqtrue_eps[,jtau] - p$Resqtrue_eps[,jtau-1])) / 
          (VecTau[jtau] - VecTau[jtau-1]) * (VecTaue - VecTau[jtau-1]) + 
          t(MatAGE[,l]) %*% p$Resqtrue_eps[,jtau-1] ) * (VecTaue > VecTau[jtau-1]) * (VecTaue <= VecTau[jtau])
    }
    #Last quantile.
    xeps[l,] = xeps[l,] + ( t(MatAGE[,l]) %*% p$Resqtrue_eps[,p$Ntau] )*( VecTaue > VecTau[p$Ntau] )

    # Laplace tails
    xeps[l,] = xeps[l,] + ( (1 / p$b1true_eps * log(VecTaue / VecTau[1])) * (VecTaue <= VecTau[1]) ) - 
      ( 1 / p$bLtrue_eps * log((1-VecTaue) / (1-VecTau[p$Ntau])) * (VecTaue > VecTau[p$Ntau]) )
  }   
  return(xeps)
}


comp.income <- function(iniage, p){

	p$age_min = iniage
	p$age_max = iniage+10
	
	etaeps <- with(p,{

		eta = comp.eta.prob(p)
		save_eta_name <- paste('eta',age_min,'.dat',sep='')	
		with( eta, save(ieta, xeta, etaprob, etacontot, etauntot, file=save_eta_name) )
		load(save_eta_name)

		epsList  = matrix(0, nrow=nage, ncol=nsim)
		etaList  = epsList
		randeta  = epsList
		enode    = epsList
	 	
	 	set.seed(77)
	  # Construct grid for draw
	  epsdraws = comp.eps(p, nsim)
		etasim <- (1:nsim) / (1+nsim)

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
        }else{
	        enode[t,i] <- which( randeta[t,i] <= etacontot[t,enode[t-1,i],] )[1]
        }

		    etaList[t,i] <- xeta[t,enode[t,i]]
		  } 
		} 


		model= list(
		epsList    = epsList,
		etaList    = etaList)   
	}) 

	savename <- paste('cohort',p$age_min,'.dat',sep='')
	save(etaeps,p,file=savename)  

  return(etaeps)
}

runMPI <- function(p){

	ai = p$firstiniage:p$lastiniage

	require(snow)  
	cl <- makeCluster(length(ai),type='MPI')
	chainN = length(cl) 
	cat('Number of Chains: ',chainN,'\n')

	start_time = proc.time()[3] 
	vals <- parLapply(cl,ai,comp.income,p)
	cat(paste('\ntotal seconds to compute income: ' , proc.time()[3] -  start_time ))

	stopCluster(cl)
	return(vals)
}