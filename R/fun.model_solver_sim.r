source('inc.modelsolver.r')

comp.solveModel <- function(p) {
	res <- with(p,{  
		
		#Setting up shock values (discrete approximation to log normal)   
		XP = qlnorm( (0:nP)/nP,-1/2*(sigP)^2,sigP )
		XP[nP+1]=100

		XT = qlnorm( (0:nT)/nT,-1/2*(sigT)^2,sigT )
		XT[nT+1]=100

		PermVecProb = array(1/nP,nP)
		ThetaVecProb = array(1/nT,nT)

		PermVec  = array(0,nP)
		ThetaVec  = array(0,nT)
		for (i in 1:nP){
			PermVec[i] = integrate(F,XP[i],XP[i+1],sigP)$value / PermVecProb[i]
		}
		for (i in 1:nT){
			ThetaVec[i] = integrate(F,XT[i],XT[i+1],sigT)$value / ThetaVecProb[i]
		}

		ThetaMat = matrix( ThetaVec,nrow=40,ncol=nT,byrow=TRUE)
		ThetaMat = rbind( ThetaMat, matrix(1,nrow=25,ncol=nT) )

		PermMat = matrix( PermVec,nrow=40,ncol=nP,byrow=TRUE)
		PermMat = rbind( PermMat, matrix(1,nrow=25,ncol=nP) )

    AlphaVec=exp(exp(exp(seq(log(log(log(aMin+1)+1)+1),log(log(log(aMax+1)+1)+1),l=n))-1)-1)-1  # set up triple exponential grid
		AlphaVec = c(AlphaVec,aHuge)

		## Profile of income level 
		Incomec = c(1.0304073e+001,1.0343448e+001,1.0382479e+001,1.0421080e+001,1.0459171e+001,1.0496674e+001,1.0533513e+001,1.0569617e+001,1.0604917e+001,1.0639347e+001,1.0672845e+001,1.0705351e+001,1.0736808e+001,1.0767163e+001,1.0796365e+001,1.0824368e+001,1.0851126e+001,1.0876600e+001,1.0900750e+001,1.0923542e+001,1.0944944e+001,1.0964928e+001,1.0983466e+001,1.1000538e+001,1.1016122e+001,1.1030203e+001,1.1042768e+001,1.1053805e+001,1.1063308e+001,1.1071273e+001,1.1077697e+001,1.1082585e+001,1.1085939e+001,1.1087769e+001,1.1088086e+001,1.1086903e+001,1.1084239e+001,1.1080114e+001,1.1074551e+001,1.1067577e+001,1.1059222e+001,1.0369837e+001,1.0369227e+001,1.0368617e+001,1.0368007e+001,1.0367397e+001,1.0366786e+001,1.0366176e+001,1.0365566e+001,1.0364956e+001,1.0364345e+001,1.0363735e+001,1.0363125e+001,1.0362515e+001,1.0361904e+001,1.0361294e+001,1.0360684e+001,1.0360074e+001,1.0359463e+001,1.0358853e+001,1.0358243e+001,1.0357633e+001,1.0357023e+001,1.0356412e+001,1.0355802e+001,1.0355192e+001)
		IncomeLevelc = exp(Incomec)

		# Income growth factors
		GList =  IncomeLevelc[-1]/IncomeLevelc[-length(IncomeLevelc)]
	

		# Construct matrix of interpolation data
		C = 0:(n+1) 
		M = 0:(n+1) 

		for l in (1:PeriodsToSolve){
			Beta      = (GList(length(GList)-l+1)^(1 - rho)*beta
			R         = Rhat/GList(length(GList)-l+1);
			ThetaVec  = ThetaMat[65-l+1,]
			PermVec   = PermMat[65-l+1,]
			
			# calculate ct from each grid point in AlphaVec
			ChiVec=nP(GothVP(AlphaVec),Rho) # inverse Euler equation
 			MuVec  = AlphaVec+ChiVec
  		M      = cbind(M, c(0,MuVec))                  # Matrix of interpolation data
  		C      = cbind(C, c(0,ChiVec))                 # Matrix of interpolation data


		}

	}) 

  return(res)
}