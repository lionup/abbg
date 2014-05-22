source('inc.modelsolver_nl_st.r')

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