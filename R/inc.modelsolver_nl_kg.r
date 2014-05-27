trans.matrix <- function(x1, x2, prob=T){
    tt <- table( x1, x2 )
    if(prob) tt <- tt / rowSums(tt)
    tt
}

econdCDF <- function(ne, wprob) {
  QcondCDF <- wprob
  for (i in 1:ne){
    for (j in 2:ne){
      QcondCDF[i,j] <- QcondCDF[i,j-1] + QcondCDF[i,j]
    }
  }
  # Force QcondCDF[,p$nw] = 1 numerically
  #QcondCDF <- QcondCDF / QcondCDF[,p$nw]   
  return(QcondCDF)
}

comp.eta.sim <- function(p){

  aa_ref = p$age_min
  #V_draw <- runif(N)  
  Vgrid<- (1:p$N) / (1+p$N)
  V_draw <- array(0, dim = c(p$nage,p$N))

  for (i in 1:p$nage){
    V_draw[i,] = sample(Vgrid)
  }

  MatAGE1 <- rep(0, p$K3+1)
  for (kk3 in 0:p$K3){    
    MatAGE1[kk3+1]=hermite( (aa_ref-p$meanAGE)/p$stdAGE, kk3 )
  }
  Mateta_true = array(0, dim = c(p$nage,p$N))

  Mateta_true[1,] = (MatAGE1 %*% p$Resqtrue_e0[,1]) * (V_draw[1,] <= p$Vectau[1])
  for (jtau in 2:p$Ntau){
    Mateta_true[1,] = Mateta_true[1,] + ( (MatAGE1 %*% 
      (p$Resqtrue_e0[,jtau]-p$Resqtrue_e0[,jtau-1]))/(p$Vectau[jtau]-p$Vectau[jtau-1]) *
      (V_draw[1,]-p$Vectau[jtau-1]) + MatAGE1 %*% p$Resqtrue_e0[,jtau-1] ) * 
      (V_draw[1,]>p$Vectau[jtau-1]) * (V_draw[1,]<=p$Vectau[jtau])
  }
  Mateta_true[1,] = Mateta_true[1,] +
    (MatAGE1 %*% p$Resqtrue_e0[,p$Ntau])*(V_draw[1,] > p$Vectau[p$Ntau])

  Mateta_true[1,]=Mateta_true[1,]+( (1/p$b1true_e0*log(V_draw[1,]/p$Vectau[1]))*(V_draw[1,]<=p$Vectau[1]) - 
    (1/p$bLtrue_e0*log((1-V_draw[1,])/(1-p$Vectau[p$Ntau]))) * (V_draw[1,]>p$Vectau[p$Ntau]))

  for (jj in 1:p$nage){  #periods before retirement
      
    aa=aa_ref+(jj-1)*2 #age last period
     
    # Eta       
    if (jj <= p$nage-1){ #already have eta at period 1, need eta for the rest nage-1 period
                       #this period is jj+1
        
      Mat = array( 0, dim=c(p$N, (p$K1+1)*(p$K2+1)) )
      for (kk1 in 0:p$K1){
        for (kk2 in 0:p$K2) {
          Mat[,kk1*(p$K2+1)+kk2+1] = hermite( (Mateta_true[jj,]-p$meanY)/p$stdY, kk1 ) * 
            hermite( ((aa+2)-p$meanAGE)/p$stdAGE, kk2 )
        }
      }
        
      #V_draw <- runif(N) 
      #First quantile
      Mateta_true[jj+1,] = ( Mat %*% p$Resqtrue[,1] ) * ( V_draw[jj+1,] <= p$Vectau[1] )

      for (jtau in 2:p$Ntau){
        Mateta_true[jj+1,] = Mateta_true[jj+1,] + 
          ( (Mat %*% (p$Resqtrue[,jtau]-p$Resqtrue[,jtau-1])) / 
          (p$Vectau[jtau]-p$Vectau[jtau-1]) *
          (V_draw[jj+1,] - p$Vectau[jtau-1]) + Mat %*% p$Resqtrue[,jtau-1] ) * 
          (V_draw[jj+1,] > p$Vectau[jtau-1]) * (V_draw[jj+1,] <= p$Vectau[jtau])
      }
      Mateta_true[jj+1,]=Mateta_true[jj+1,] + (Mat %*% p$Resqtrue[,p$Ntau]) * (V_draw[jj+1,]>p$Vectau[p$Ntau])

      Mateta_true[jj+1,]=Mateta_true[jj+1,] + ( (1 / p$b1true * log(V_draw[jj+1,]/p$Vectau[1])) * (V_draw[jj+1,] <= p$Vectau[1]) - 
          (1 / p$bLtrue * log((1-V_draw[jj+1,])/(1-p$Vectau[p$Ntau]))) * (V_draw[jj+1,] > p$Vectau[p$Ntau]) )
    }
  }

  save_Mateta_name <- paste('Mateta',aa_ref,'.dat',sep='')
  save(Mateta_true, file=save_Mateta_name)

  return(Mateta_true)
}

comp.eta.prob <- function(p){
  
  # get the simulations of workers
  #Mateta_true <- comp.eta.sim(p)
  save_Mateta_name <- paste('Mateta',p$age_min,'.dat',sep='')
  load(save_Mateta_name)

  # Quantiles of eta and epsilon, by age
  xeta <- array( 0, dim=c(p$nage, p$nbin) )  #bins
  veta <- seq(1/(2*p$nbin), (2*p$nbin-1)/(2*p$nbin), l=(2*p$nbin-1)) #bins and its median
  oddnode <- seq(1,(2*p$nbin-1),2)  #income node is the median in each interval, odd node
  age = seq(p$age_min, p$age_re-2, 2)

  for (i in 1:p$nage){ #periods before retirement
    xeta[i,] <- quantile( Mateta_true[i,], veta, names=F, na.rm = T )[oddnode]
  }  

  long <- data.table( pid = 1:p$N, age=rep(age, each=p$N), income = c(t(Mateta_true)) )
  setkey(long, pid, age)
  long[, q:=as.numeric(cut_number(income, n = p$nbin)), age] #give a bin # for each person each age
  long$income <- NULL
  wide <- reshape(long, idvar='pid', timevar='age', direction='wide') #each person has all age in a row

  etaprob <- array( 0,dim=c(p$nage, p$nbin, p$nbin) ) #last period to this period
  etaprob[1,,] <- 1/p$nbin  #first period, same prob

  for (t in 2:p$nage){ #today
    x1 <- paste('q',age[t-1],sep='.')
    x2 <- paste('q',age[t],sep='.')
    etaprob[t,,] <- trans.matrix( wide[[x1]], wide[[x2]] ) #today
  }     

  etacontot <- etaprob
  for ( i in 1:p$nage ) {
    etacontot[i,,] <- econdCDF( p$nbin, etaprob[i,,] )
  }  

  etauntot <- etaprob[1,1,]
  for(i in 2:p$nbin){
    etauntot[i] <- etauntot[i-1]  + etauntot[i]  
  }

  mineta <- array( 0, dim=c(p$nage, p$nbin) ) #this period and last eta
  maxeta <- array( p$nbin+1, dim=c(p$nage, p$nbin) )
  for(t in 2:p$nage){ #today
    for(w in 1:p$nbin){   #for yesterday
      for(ww in 1:p$nbin){ #for today
        if (etaprob[t,w,ww] > 0) {
          break
        }else{
          mineta[t,w] = ww
        }
      }   

      for(ww in p$nbin:1){ #for today
        if (etaprob[t,w,ww] > 0) {
          break
        }else{
          maxeta[t,w] = ww
        }
      }   
    } 
  }    

  eta = list(
  xeta      = xeta,
  etaprob   = etaprob,
  mineta    = mineta,
  maxeta    = maxeta, 
  etacontot = etacontot,
  etauntot  = etauntot)
  
  return(eta)
}

comp.eps <- function(p, ntra){ 
  xeps <- array( 0, dim=c(p$nage, ntra) ) #income grid
  VecTaue <- seq(1/(1+p$neps), p$neps/(1+p$neps), l=ntra)
  #VecTaue <- (1:ntra) / (1+ntra) #get the quantile of interpolation node

  # generate the age vector 
  xt  <- seq( p$age_min, (p$age_min + (p$nage-1) * 2 ), by=2 )
  xtL <- ( xt - p$meanAGE ) / p$stdAGE #standardize AGE

  # Create MatAGE which is the hermite of AGE
  MatAGE <- array( 0,dim=c(p$K4+1, p$nage) )
  for (kk4 in 0:p$K4){    
    MatAGE[kk4+1,]=hermite(xtL,kk4)
  }

  for (l in 1:p$nage){
    xeps[l,] = ( t(MatAGE[,l]) %*% p$Resqtrue_eps[,1] )*( VecTaue <= p$Vectau[1] )
    for (jtau in 2:p$Ntau){
        xeps[l,] = xeps[l,] + ( (t(MatAGE[,l]) %*% (p$Resqtrue_eps[,jtau] - p$Resqtrue_eps[,jtau-1])) / 
          (p$Vectau[jtau] - p$Vectau[jtau-1]) * (VecTaue - p$Vectau[jtau-1]) + 
          t(MatAGE[,l]) %*% p$Resqtrue_eps[,jtau-1] ) * (VecTaue > p$Vectau[jtau-1]) * (VecTaue <= p$Vectau[jtau])
    }
    #Last quantile.
    xeps[l,] = xeps[l,] + ( t(MatAGE[,l]) %*% p$Resqtrue_eps[,p$Ntau] )*( VecTaue > p$Vectau[p$Ntau] )

    # Laplace tails
    xeps[l,] = xeps[l,] + ( (1 / p$b1true_eps * log(VecTaue / p$Vectau[1])) * (VecTaue <= p$Vectau[1]) ) - 
      ( 1 / p$bLtrue_eps * log((1-VecTaue) / (1-p$Vectau[p$Ntau])) * (VecTaue > p$Vectau[p$Ntau]) )
  }   
  return(xeps)
}


#CRRA marginal utility function      
uP <- function(c,Rho){
  return(c^(-Rho))
}

# Gothic V prime function
GothVP <- function(p, a, inc, FC, FM){
  mtp = p$R*a + inc  # money next period
  EUP = uP( Cnextp(mtp,FC,FM),p$rho )    
  return(EUP)
}

Cnextp <- function(m,ctp1,mtp1){
# Cnextp is constructed by interpolation to be the next-period consumption function Ct+1()
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


#####################################################
FindLinProb1 <- function(xi,x){
  #this takes in xi, and finds two points in either side of it in x
  #and returns the indices of them and associated probabilities 
  nx = length(x)
  LocL = tail(which(x<=xi),1)
  fmin = c(1,1,0) #xi<=x[1]

  if(xi>=x[nx]){
    fmin[1] = nx-1
    fmin[2] = 0
    fmin[3] = 1
  }else if(xi>x[1]){
    fmin[1] = LocL 
    inteta <- (2*xi - x[LocL+1]-x[LocL]) / (x[LocL+1]-x[LocL])
    fmin[2] = 0.5 * (1 - inteta) 
    fmin[3] = 1 - fmin[2]
  }
  return(fmin)
}