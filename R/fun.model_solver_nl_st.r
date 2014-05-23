source('inc.modelsolver_nl_st.r')

comp.income <- function(iniage, p){

	require(EQL)
	require(data.table)
	require(ggplot2)

	p$age_min = iniage
	p$age_max = iniage+10

	#eta = comp.eta.prob(p)
	save_eta_name <- paste('eta',p$age_min,'.dat',sep='')	
	#with( eta, save(ieta, xeta, etaprob, etacontot, etauntot, file=save_eta_name) )
	load(save_eta_name)

	epsList  = matrix(0, nrow=p$nage, ncol=p$nsim)
	etaList  = epsList
	randeta  = epsList
	enode    = epsList
 	
 	set.seed(77)
  # Construct grid for draw

  cat('\n after load',iniage,file="help.txt",append=T)
  epsdraws = comp.eps(p, p$nsim)
  
 	cat('\n before loop',iniage,file="help.txt",append=T)

	etasim <- (1:p$nsim) / (1+p$nsim)

	#loop over life cycle
	for (t in 1:p$nage) {  
		# Sample randomly from eps grid to get current eps
		epsList[t,] = sample(epsdraws[t,]) 

		# generate random draw on unit interval for current eta
		randeta[t,] = sample(etasim)

		cat('\n before inner loop time',t,'age',iniage,file="help.txt",append=T)

		#loop for different individuals
		for (i in 1:p$nsim){
			# find persistent income draw
	  	if(t==1){  
        enode[1,i] <- which( randeta[t,i] <= etauntot )[1]
      }else{
        enode[t,i] <- which( randeta[t,i] <= etacontot[t,enode[t-1,i],] )[1]
      }

	    etaList[t,i] <- xeta[t,enode[t,i]]
	  } 
	} 

	cat('\n after loop',iniage,file="help.txt",append=T)

	model= list(epsList = epsList, etaList = etaList)   

	savename <- paste('cohort',p$age_min,'.dat',sep='')
	save(model,p,file=savename)  

  return(iniage)
}

runMPI <- function(p){

	ai = p$firstiniage:p$lastiniage

	if (p$mode=='mpi') {

		cat('[mode=mpi] USING MPI !!!!! \n')
		require(snow)  
		cl <- makeCluster(length(ai),type='MPI')
		chainN = length(cl) 
		cat('Number of Chains: ',chainN,'\n')

			clusterExport(cl, list("trans.matrix", 
			"econdCDF",
			"comp.eta.sim", 
			"comp.eta.prob", 
			"comp.eps",
			"hermite"))

		start_time = proc.time()[3] 
		vals <- parSapply(cl,ai,comp.income,p)
		cat(paste('\ntotal seconds to compute cohort: ' , proc.time()[3] -  start_time ))
		stopCluster(cl)

	} else if (p$mode == 'serial') {
    cat('[mode=serial] NOT USING MPI !!!!! \n')
    start_time = proc.time()[3]   
    vals <- sapply(ai,comp.income,p)
    cat(paste('\ntotal seconds to compute cohort: ' , proc.time()[3] -  start_time ))
  }  

	return(vals)
}