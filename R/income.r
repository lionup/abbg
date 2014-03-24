rm(list=ls())
#require(R.matlab)
require(EQL)
#require(data.table)
setwd('~/git/abbg/R')
#data <- readMat('data_hermite.mat')
load('mat.dat')
source('inc.modelsolver.r')

# two periods
#declare parameters
p <- list()
p$ne <- 23  #eta permanent income
p$ntau <- 11  #number of tau used as interpolation node
p$nt <- 6  #from 30 ~ 65
p$K1  <- 3  #third order hermite for income
p$K2  <- 1  #first order hermite for age
p$nsim  = 10000

xe <- array(0, dim=c(p$nt, p$ne)) #income grid
VecTau <- (1:p$ntau)/(1+p$ntau) #get the quantile of interpolation node
VecTaue <- (1:p$ne)/(1+p$ne) #get the quantile of interpolation node

# generate the age vector 
xt  <- 30:(30+p$nt-1)
xtL <- ( xt - data$meanAGE ) / data$stdAGE #standardize AGE
#xtL <- ( xt - mean(xt) ) / sd(xt) #standardize AGE

# Create MatAGE which is the hermite of AGE
MatAGE <- array( 0,dim=c(p$K2+1, p$nt) )
for (kk2 in 0:p$K2){    
  MatAGE[kk2+1,]=hermite(xtL,kk2)
}

# get the income bounds from the data
#ageY <- data.table(age=as.numeric(data$AGE), income=as.numeric(data$Y))
#setkey(ageY, age)
#minY <- rep( 0, p$nt)
#maxY <- rep( 0, p$nt)
#for (i in xt){
#  minY[i-29] = min(ageY[J(i)]$income)
#  maxY[i-29] = max(ageY[J(i)]$income)
#}

#  create the grid nodes of eta_1
# Resqfinal.e0 is the a0 and a1 for age1. 
#xe[1,] <- t(MatAGE[,1]) %*% data$Resqfinal.e0
#if ( xe[1,1]   > minY[1] ) minY[1] = xe[1,1]
#if ( xe[1,p$ne]< maxY[1] ) maxY[1] = xe[1,p$ne]
xe[1,] = ( t(MatAGE[,1]) %*% data$Resqfinal.e0[,1] )*( VecTaue <= VecTau[1] )
for (jtau in 2:p$ntau){
    xe[1,] = xe[1,] + ( (t(MatAGE[,1]) %*% (data$Resqfinal.e0[,jtau] - data$Resqfinal.e0[,jtau-1])) / 
    	(VecTau[jtau] - VecTau[jtau-1]) * (VecTaue - VecTau[jtau-1]) + 
    	t(MatAGE[,1]) %*% data$Resqfinal.e0[,jtau-1] ) * (VecTaue > VecTau[jtau-1]) * (VecTaue <= VecTau[jtau])
}
#Last quantile.
xe[1,] = xe[1,] + ( t(MatAGE[,1]) %*% data$Resqfinal.e0[,p$ntau] )*( VecTaue > VecTau[p$ntau] )

# Laplace tails
xe[1,] = xe[1,] + ( (1 / data$b1.e0 * log(VecTaue / VecTau[1])) * (VecTaue <= VecTau[1]) ) - 
	( 1 / data$bL.e0 * log((1-VecTaue) / (1-VecTau[p$ntau])) * (VecTaue > VecTau[p$ntau]) )

# create the transition matrix 
# row is this 
etaprob <- array( 0,dim=c(p$nt, p$ne, p$ne) )
# the prob of first period is uniform
etaprob[1,,] <- 1/p$ne

#transform eta_(t-1) into hermite
for (t in 2:p$nt){
  xeL <- ( xe[t-1,] - data$meanY ) / data$stdY #standardize eta_1
  #xeL <- ( xe[t-1,] - mean(xe[t-1,]) ) / sd(xe[t-1,]) 
  #xeL <- xe[t-1,]
	
	#evaluate eta_(t-1) by third order hermite
	Mat <- array( 0, dim=c( (p$K1+1)*(p$K2+1), p$ne) )
	for (kk1 in 0:p$K1){   
		for (kk2 in 0:p$K2){  
	    Mat[kk1*(p$K2+1)+kk2+1,] = hermite(xeL,kk1) * MatAGE[kk2+1,t]
		}
	}

	#for each grid node of eta_(t-1), locate the quantiles for all grid nodes of eta_t
	#for each grid node of eta_(t-1), find its value in period t on each tau from vectau
	#multiply the hermite of eta_(t-1) by ak, ne*ne matrix, row is e this period, col is e next period 
	#Resqfinal is ak
	xee <- t(Mat) %*% data$Resqfinal

	#  get the min and mean nodes on the grid of eta_2
  # min is the extrapolate of vectau[1] node of xe[t-1,1]
  # max is the extrapolate of vectau[p$ne] node of xe[t-1,p$ne]
  #mine <- max(-6, min(xee))
  #meane <- mean(xee)  
  #flip the mine w.r.t meane and get the maxe 
  #maxe <- meane - mine
  #maxe <- min(6, max(xee))
  #mine = ( t(Mat[,1]   ) %*% data$Resqfinal[,1] )
  #maxe = ( t(Mat[,p$ne]) %*% data$Resqfinal[,p$ntau] )
  
  # Laplace tails
  mine = xee[1]           + 1 / data$b1 * log( VecTaue[1]        / VecTau[1]          )  
  mine = max(-6, mine)
  maxe = xee[length(xee)] - 1 / data$bL * log( (1-VecTaue[p$ne]) / (1-VecTau[p$ntau]) ) 
  maxe = min( 6, maxe)
  
  
	# then discretize the grid into ne points
	xe[t,] <- seq(mine, maxe, l=p$ne)
  
  #find the probabilities for transitions from lw(t-1) to lw(t)
  #for each grid node of eta_(t-1), find the tau of all grid note of eta_t
  tauee <- array( 0, dim=c(p$ne, p$ne) )
	for (i in 1:p$ne){ #eta_(t-1)
		for (j in 1:p$ne) { #eta_t
			#if eta_t is in the range of xee[i,], interpolate
			if( xe[t,j] <= xee[i,p$ntau] & xe[t,j] >= xee[i,1] ){
        tauee[i,j] <- approx(xee[i,],VecTau, xe[t,j])$y
			}else if( xe[t,j] < xee[i,1]){ #if eta_t is out of  the range of xee[i,], find tau using opt
			  tauee[i,j] <-      VecTau[1]         * exp( data$b1 * (xe[t,j]       - xee[i,1]) )
			}else{
			  tauee[i,j] <- 1 - (1-VecTau[p$ntau]) * exp( data$bL * (xee[i,p$ntau] - xe[t,j]) )
			}      
		}		
	}
  #cutting the unit interval into 11 portions by the quantiles calculated 
  #the transition probability is the length of each protions
  # a, b,c => a+(b-a)/2; (b-a)/2+(c-b)/2; 1-c+(c-b)/2
  tau_1 = tauee[,-p$ne] #23*22 1:22
  tau_2 = tauee[,-1]  #23*22 2:23
  tau_3 = (tau_2-tau_1)/2 #2-1 ~23-22
  etaprob[t,,1] = tauee[,1] + tau_3[,1]
  etaprob[t,,2:(p$ne-1)] = tau_3[,-(p$ne-1)] + tau_3[,-1]
  etaprob[t,,p$ne]  = 1 - tauee[,p$ne] + tau_3[,(p$ne-1)]
 
}
etacontot <- etaprob
for (i in 1:p$nt) {
  etacontot[i,,] <- econdCDF(p, etaprob[i,,])
}  
etauntot <- etaprob[1,1,]
for(i in 2:p$ne){
	etauntot[i] <- etauntot[i-1]  + etauntot[i]  
}

#simulate permanent income
etalong <- rep(0, p$nt)
etaind <- array(0, dim=c(p$nsim, p$nt))
randeta <- array(0, dim=c(p$nsim, p$nt))
visite  <- array(0, dim=c(p$nt, p$ne))
set.seed(123)
for (i in 1:p$nt){
	randeta[,i] <- runif(p$nsim)  #income
}

#loop for different individuals
for (i in 1:p$nsim){
 
	#loop over life cycle
	for (t in 1:p$nt) {     

 		# find income draw for individual i at period 1
		if(t==1){
			enode <- which( randeta[i,t] <= etauntot )[1]
			visite[t,enode] = visite[t,enode] + 1
		# find income draw for individual i at period t
		}else{
			enode <- which( randeta[i,t] <= etacontot[t,enode,] )[1]
			visite[t,enode] = visite[t,enode] + 1
		}

		etaind[i,t]<-xe[t,enode]
		etalong[t] <- etalong[t] + xe[t,enode]
	}		
}	