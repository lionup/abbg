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

    if(age_min==30){ 
      save(Mateta_true, file='Mateta_even.dat')
    }else{
      save(Mateta_true, file='Mateta_odd.dat')
    }  
    Mateta_true
  }) 
  return(res)
}

comp.eta.prob <- function(p){
  res <- with(p,{ 

    # get the simulations of workers
    #Mateta_true <- comp.eta.sim(p)
    if(age_min==30){ 
      load('Mateta_even.dat')
    }else{
      load('Mateta_odd.dat')
    }  

    # Quantiles of eta and epsilon, by age
    xeta <- array( 0, dim=c(nage+ntr, nbin) )  #bins
    veta <- seq(1/(2*nbin), (2*nbin-1)/(2*nbin), l=(2*nbin-1)) #bins and its median
    oddnode <- seq(1,(2*nbin-1),2)  #income node is the median in each interval, odd node
    age = seq(age_min, age_re-2, 2)

    for (i in 1:nage){ #periods before retirement
      xeta[i,] <- quantile( Mateta_true[i,], veta, names=F, na.rm = T )[oddnode]
    }  

    for( i in (nage+1):(nage+ntr) ){ #periods after retirement
      xeta[i,] <- xeta[nage,]
    }  

    long <- data.table( pid = 1:N, age=rep(age, each=N), income = c(t(Mateta_true)) )
    setkey(long, pid, age)
    long[, q:=as.numeric(cut_number(income, n = nbin)), age] #give a bin # for each person each age
    long$income <- NULL
    wide <- reshape(long, idvar='pid', timevar='age', direction='wide') #each person has all age in a row

    etaprob <- array( 0,dim=c(nage+ntr, nbin, nbin) ) #last period to this period
    etaprob[1,,] <- 1/nbin  #first period, same prob

    for (t in 2:nage){ #today
      x1 <- paste('q',age[t-1],sep='.')
      x2 <- paste('q',age[t],sep='.')
      etaprob[t,,] <- trans.matrix( wide[[x1]], wide[[x2]] ) #today
    }  

    for( t in (nage+1):(nage+ntr) ){
      etaprob[t,,] <- diag(nbin)
    }  

    mineta <- array( 0, dim=c(nage+ntr, nbin) ) 
    maxeta <- array( nbin+1, dim=c(nage+ntr, nbin) )
    for ( t in 2:(nage+ntr) ){ #today
      for(w in 1:nbin){   #for yesterday
        for(ww in 1:nbin){ #for today
          if (etaprob[t,w,ww] > 0) {
            break
          }else{
            mineta[t,w] = ww
          }
        }   

        for(ww in nbin:1){ #for today
          if (etaprob[t,w,ww] > 0) {
            break
          }else{
            maxeta[t,w] = ww
          }
        }   
      } 
    }    

    etacontot <- etaprob
    for ( i in 1:(nage+ntr) ) {
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
    mineta    = mineta,
    maxeta    = maxeta, 
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

#CRRA marginal utility function      
uP <- function(c,Rho){
  return(c^(-Rho))
}

# Gothic V prime function
GothVP <- function(p, t, a, inc, theta, perm, thetaP, permP, FC, FM, minperm, maxperm){
  EUP = rep( 0, length(a) ) 
  for ( j in (minperm+1):(maxperm-1) ){
    if (t>p$nage){  #retired no income shock   
      mtp = p$R*a + inc * perm[j]  # money next period
      EUP = EUP + uP( Cnextp(mtp,FC[j,],FM[j,]),p$rho )* permP[j]       
    }else{
      for ( i in 1:p$neps ){
        mtp = p$R*a + inc * perm[j] * theta[i]  # money next period
        EUP = EUP + uP( Cnextp(mtp,FC[j,],FM[j,]),p$rho )*thetaP[i]*permP[j]       
      }
    }
  }  
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

