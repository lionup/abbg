
# integrate x  by log normal
F  <- function(x, Sigma) {
	x * dlnorm(x,-1/2*(Sigma)^2,Sigma)
}	

DiscreteApproxToMeanOneLogNormal <- function(std,numofshockpoints){
  LevelAdjustingParameter = -(1/2)*(std)^2
  ListOfEdgePoints = qlnorm( (0:numofshockpoints)/numofshockpoints,LevelAdjustingParameter,std )
  ListOfEdgePoints[numofshockpoints+1]=100

  shocklist  = rep(0,numofshockpoints)
  for (i in 1:numofshockpoints){
    shocklist[i] = integrate(F,ListOfEdgePoints[i],ListOfEdgePoints[i+1],std)$value * numofshockpoints
  }
  return(shocklist)
}

#CRRA marginal utility function      
uP <- function(c,Rho){
  return(c^(-Rho))
}

# Gothic V prime function
GothVP <- function(p, a, G, theta, perm, thetaP, permP, C, M){
  EUP = rep( 0, length(a) ) 
  for ( i in 1:length(theta) ){
    for (j in 1:length(perm) ){
      mtp = p$R*a/(G*perm[j]) + theta[i]  # money next period
      EUP = EUP + uP( perm[j]*Cnextp(mtp,C,M),p$rho )*thetaP[i]*permP[j]       
    }
  }
  EUP = EUP * p$beta * p$R * G^(-p$rho) 
  return(EUP)
}

# Ct+1 function 
Cnextp <- function(m,C,M){
# Cnextp is constructed by interpolation to be the next-period consumption function Ct+1()
  mtp1 = M[,ncol(M)]  # data for the next-period consumption function
  ctp1 = C[,ncol(C)]  # data for the next-period consumption function

  c = rep( 0, length(m) ) 
  end = length(mtp1)

  # extrapolate above maximal m
  iAbove = m >= mtp1[end]
  slopeAbove  = ( ctp1[end]-ctp1[end-1] ) / ( mtp1[end]-mtp1[end-1] )
  c[iAbove]   = ctp1[end] + ( m[iAbove]-mtp1[end] )*slopeAbove

  # extrapolate below minimal m
  iBelow = m <= mtp1[1]
  slopeBelow  = 1
  c[iBelow]   = ctp1[1] + ( m[iBelow]-mtp1[1] )*slopeBelow

  # interpolate
  iInterp = !(iAbove | iBelow)
  c[iInterp]  = approx( mtp1,ctp1,m[iInterp] )$y
  return(c)  
}

