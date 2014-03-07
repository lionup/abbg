source('inc.modelsolver.r')

comp.solveModel <- function(p) {
	res <- with(p,{  
		
		#Setting up shock values (discrete approximation to log normal)   
		X = qlnorm((0:NumOfThetaShockPoints)/NumOfThetaShockPoints,-1/2*(Sigma)^2,Sigma)
		X[NumOfThetaShockPoints+1]=1e3

		ThetaProb = array(1/NumOfThetaShockPoints,NumOfThetaShockPoints)
		ThetaVals = array(0,NumOfThetaShockPoints)
		for (i in 1:NumOfThetaShockPoints){
			ThetaVals[i] = integrate(F,X[i],X[i+1],Sigma)$value/ThetaProb[i]
		}




	}) 

  return(res)
}