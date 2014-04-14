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

comp.eta.sim <- function(p, N, age){
  res <- with(p,{  
    # Number of individuals
    #N=999999
    #N=1000
    aa_ref=age
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

    for (jj in 1:nage){
        
      aa=aa_ref+(jj-1)*2
       
      # Eta       
      if (jj <= nage-1){ #this period is jj+1
          
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

    save(Mateta_true, file='Mateta.dat')
    Mateta_true
  }) 
  return(res)
}

comp.eta.prob <- function(p, N){
  res <- with(p,{ 

    # get the simulations of workers
    #Mateta_true <- comp.eta.sim(p,N)
    load('Mateta.dat')

    # Quantiles of eta and epsilon, by age
    xeta <- array( 0, dim=c(nage, nbin) )  #23 bin
    veta <- seq(1/(2*nbin), (2*nbin-1)/(2*nbin), l=(2*nbin-1)) #45 nodes
    oddnode <- seq(1,(2*nbin-1),2)  #income node is the even node
    age = seq(30, (30+(nage-1)*2), 2)

    for (i in 1:nage){
      xeta[i,] <- quantile( Mateta_true[i,], veta, names=F, na.rm = T )[oddnode]
    }    

    long <- data.table( pid = 1:N, age=rep(age, each=N), income = c(t(Mateta_true)) )
    setkey(long, pid, age)
    long[, q:=as.numeric(cut_number(income, n = nbin)), age]
    long$income <- NULL
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
  xt  <- seq( 30, (30 + (p$nage-1) * 2 ), by=2 )
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
GothVP <- function(p, a, theta, perm, thetaP, permP, FC, FM){
  EUP = rep( 0, length(a) ) 
  for ( i in 1:p$neps ){
    for ( j in 1:p$nbin ){
      mtp = p$R*a + perm[j] * theta[i]  # money next period
      EUP = EUP + uP( Cnextp(mtp,FC[j,],FM[j,]),p$rho )*thetaP[i]*permP[j]       
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

