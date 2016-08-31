trans.matrix <- function(x1, x2, prob=T){
    tt <- table( x1, x2 )
    if(prob) tt <- tt / rowSums(tt)
    tt
}

#####################################################
comp.eta.sim <- function(p){

  Vgrid <- (1:p$N) / (1+p$N)
  V_draw <- array(0, dim = c(p$Twork,p$N))
  for (i in 1:p$Twork) {
    V_draw[i,] = sample(Vgrid)
  }

  aa_ref = p$age_min
  MatAGE1 <- rep(0, p$K3+1)
  for (kk3 in 0:p$K3){
    MatAGE1[kk3+1] <- hermite( (aa_ref-p$meanAGE)/p$stdAGE, kk3 )
  }

  Mateta <- array(0, dim = c(p$Twork,p$N))
  Mateta[1,] <- (MatAGE1 %*% p$Resqtrue_e0[,1]) * (V_draw[1,] <= p$Vectau[1])
  for (jtau in 2:p$Ntau){
    Mateta[1,] <- Mateta[1,] +
      ( (MatAGE1 %*% (p$Resqtrue_e0[,jtau] - p$Resqtrue_e0[,jtau-1])) /
      (p$Vectau[jtau] - p$Vectau[jtau-1]) *
      (V_draw[1,] - p$Vectau[jtau-1]) + MatAGE1 %*% p$Resqtrue_e0[,jtau-1] ) *
      (V_draw[1,] > p$Vectau[jtau-1]) * (V_draw[1,] <= p$Vectau[jtau])
  }
  Mateta[1,] <- Mateta[1,] + (MatAGE1 %*% p$Resqtrue_e0[,p$Ntau])*(V_draw[1,] > p$Vectau[p$Ntau])

  Mateta[1,] <- Mateta[1,]+( (1/p$b1true_e0*log(V_draw[1,]/p$Vectau[1]))*(V_draw[1,]<=p$Vectau[1]) -
    (1/p$bLtrue_e0*log((1-V_draw[1,])/(1-p$Vectau[p$Ntau]))) * (V_draw[1,]>p$Vectau[p$Ntau]))

  for (jj in 1:p$Twork){  #periods before retirement

    aa=aa_ref+(jj-1)*2 #age last period

    # Eta
    if (jj <= p$Twork-1){ #already have eta at period 1, need eta for the rest Twork-1 period
                       #this period is jj+1
      Mat <- array( 0, dim=c(p$N, (p$K1+1)*(p$K2+1)) )
      for (kk1 in 0:p$K1){
        for (kk2 in 0:p$K2) {
          Mat[,kk1*(p$K2+1)+kk2+1] <- hermite( (Mateta[jj,]-p$meanY)/p$stdY, kk1 ) *
            hermite( ((aa+2)-p$meanAGE)/p$stdAGE, kk2 )
        }
      }

      #V_draw <- runif(N)
      #First quantile
      Mateta[jj+1,] <- ( Mat %*% p$Resqtrue[,1] ) * ( V_draw[jj+1,] <= p$Vectau[1] )

      for (jtau in 2:p$Ntau){
        Mateta[jj+1,] <- Mateta[jj+1,] +
          ( (Mat %*% (p$Resqtrue[,jtau]-p$Resqtrue[,jtau-1])) /
          (p$Vectau[jtau]-p$Vectau[jtau-1]) *
          (V_draw[jj+1,] - p$Vectau[jtau-1]) + Mat %*% p$Resqtrue[,jtau-1] ) *
          (V_draw[jj+1,] > p$Vectau[jtau-1]) * (V_draw[jj+1,] <= p$Vectau[jtau])
      }
      Mateta[jj+1,] <- Mateta[jj+1,] + (Mat %*% p$Resqtrue[,p$Ntau]) * (V_draw[jj+1,]>p$Vectau[p$Ntau])

      Mateta[jj+1,] <- Mateta[jj+1,] + ( (1 / p$b1true * log(V_draw[jj+1,]/p$Vectau[1])) * (V_draw[jj+1,] <= p$Vectau[1]) -
          (1 / p$bLtrue * log((1-V_draw[jj+1,])/(1-p$Vectau[p$Ntau]))) * (V_draw[jj+1,] > p$Vectau[p$Ntau]) )
    }
  }

  save_Mateta_name <- paste('Mateta',aa_ref,'.dat',sep='')
  save(Mateta, file=save_Mateta_name)

  return(Mateta)
}

#####################################################
comp.eta.prob <- function(p){
  res <- with(p,{
    #attach(p)

    # get the simulations of workers
    #Mateta <- comp.eta.sim(p)
    save_Mateta_name <- paste('Mateta',age_min,'.dat',sep='')
    load(save_Mateta_name)

    # Quantiles of eta and epsilon, by age
    zgrid <- array( 0, dim=c(Twork, ngpz) )  #bins
    ztrans <- array( 0,dim=c(Twork-1, ngpz, ngpz) ) #last period to this period

    veta <- seq(1/(2*ngpz), (2*ngpz-1)/(2*ngpz), l=ngpz) #median of each bin
    age = 1:Twork

    for (i in 1:Twork){ #periods before retirement
      zgrid[i,] <- quantile( Mateta[i,], veta, names=F, na.rm = T )
    }

    long <- data.table( pid = 1:N, age=rep(age, each=N), income = c(t(Mateta)) )
    setkey(long, pid, age)
    long[, q:=as.numeric(cut_number(income, n = ngpz)), age] #give a bin # for each person each age
    long$income <- NULL
    wide <- reshape(long, idvar='pid', timevar='age', direction='wide') #each person has all age in a row

    for ( t in 1:(Twork-1) ){ #today
      x1 <- paste('q',age[t  ],sep='.')
      x2 <- paste('q',age[t+1],sep='.')
      ztrans[t,,] <- trans.matrix( wide[[x1]], wide[[x2]] ) #today
    }

    zdist <- 1/ngpz #unconditonal distribution, each bin is 0.01 of the population, so each node is with prob 0.01.

    #find variance
    lvar <- rep(0,Twork)
    for (it in 1:Twork) {
      lvar[it] <- sum(zgrid[it,]^2 * zdist) - sum(zgrid[it,] * zdist)^2
    }

    rr=list(zdist=zdist, zgrid=zgrid, ztrans=ztrans, varzapprox=lvar)
  })
  return(res)
}

#####################################################
comp.eps <- function(p){
  xeps <- array( 0, dim=c(p$Twork, p$ngpe) ) #income grid
  VecTaue <- seq(1/(2*p$ngpe), (2*p$ngpe-1)/(2*p$ngpe), l=p$ngpe)  #get the quantile of interpolation node

  # generate the age vector
  xt  <- seq( p$age_min, p$age_re-2, by=2 )
  xtL <- ( xt - p$meanAGE ) / p$stdAGE #standardize AGE

  # Create MatAGE which is the hermite of AGE
  MatAGE <- array( 0,dim=c(p$K4+1, p$Twork) )
  for (kk4 in 0:p$K4){
    MatAGE[kk4+1,]=hermite(xtL,kk4)
  }

  for (l in 1:p$Twork){
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


#####################################################
#Ygross - tax function - Ynet
FnGrossInc <- function(lx,lnet,p,lstax){
  #lx is gross, lnet is net
  lf = (1-p$pentax - p$btax)*lx + p$btax*( lx^(-p$ptax) + lstax )^(-1/p$ptax) - lnet
}

#####################################################
FnTaxParamNet <- function(lstax, p, kappa, popsize, zgrid, egrid, zdist, edist, moment=TRUE){
  res <- with(p,{
    if (Display==1) cat("Tax function parameter: ", lstax, "\n")

    avearnspre     <- rep( 0, Twork )
    #avearnspost    <- avearnspre
    #avlearnspre    <- avearnspre
    #avlearnspost   <- avearnspre
    #avearnspre2    <- avearnspre
    #avearnspost2   <- avearnspre
    #avlearnspre2   <- avearnspre
    #avlearnspost2  <- avearnspre

    ltottax <- 0.0
    ltotlabincpre <- 0.0
    #ltotlabincpost = 0.0

    ygrid    <- array( 0, dim=c(Twork,ngpz,ngpe) ) #earnings grid
    ypregrid <- ygrid

    for(it in 1:Twork){
      for(iz in 1:ngpz){
        for(ie in 1:ngpe){
          ygrid[it,iz,ie] = exp( kappa[it] + zgrid[it,iz] + egrid[it,ie] )

          #get implied gross income at this point, to use in constructtion of soc sec system
          lygross <- ygrid[it,iz,ie]/(1.0-pentax-btax)
          ypregrid[it,iz,ie] <- uniroot(FnGrossInc, c(0, lygross*3), ygrid[it,iz,ie], p, lstax,
            extendInt="yes", tol=1, maxiter=30)$root

          ltotlabincpre <- ltotlabincpre + ypregrid[it,iz,ie]*zdist*edist*popsize[it]
          #ltotlabincpost = ltotlabincpost + ygrid[it,iz,ie]*zdist*edist*popsize[it]
          ltottax <- ltottax + btax*( ypregrid[it,iz,ie] - (ypregrid[it,iz,ie]^(-ptax) + lstax)^(-1/ptax) ) *
            zdist*edist*popsize[it]

          avearnspre[it] <- avearnspre[it] + ypregrid[it,iz,ie]*zdist*edist
          #avearnspre2[it] = avearnspre2[it] + (ypregrid[it,iz,ie]^2)*zdist*edist
          #avearnspost[it] = avearnspost[it] + ygrid[it,iz,ie]*zdist*edist
          #avearnspost2[it] = avearnspost2[it] + (ygrid[it,iz,ie]^2)*zdist*edist
          #avlearnspre[it] = avlearnspre[it] + log(ypregrid[it,iz,ie])*zdist*edist
          #avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it,iz,ie]^2)*zdist*edist
          #avlearnspost[it] = avlearnspost[it] + log(ygrid[it,iz,ie])*zdist*edist
          #avlearnspost2[it] = avlearnspost2[it] + log(ygrid[it,iz,ie]^2)*zdist*edist
        }
      }
    }

    #varearnspre    = avearnspre2     - avearnspre^2
    #varearnspost   = avearnspost2   - avearnspost^2
    #varlearnspre   = avlearnspre2   - avlearnspre^2
    #varlearnspost  = avlearnspost2 - avlearnspost^2

    if (Display==1) cat('Tax revenue / Pre-tax labor income: ', (ltottax/ltotlabincpre)*100, '%\n')

    if(moment){
      rr <- ltottax/ltotlabincpre - targetTaxToLabinc
    }else{
      rr <- list(ygrid = ygrid,
                ypregrid = ypregrid,
                avearnspre   = avearnspre)
    }

  })
  return(res)
}

#####################################################
#tax function
FnTax <- function(ly,p){
  lf = p$btax*( ly - (ly^(-p$ptax) + p$stax)^(-1/p$ptax) ) + p$pentax*ly
}

#####################################################
FnSSParam <- function(lsspar, p, avearnspre, mgrid, moment=TRUE){
  res <- with(p,{
    if (Display==1) cat('SS benefit parameter: ', lsspar, '\n')

    lavearns <- sum(avearnspre)/Twork #cross-sectional average gross earnings
    ssbend1  <- 0.18*lavearns
    ssbend2  <- 1.1*lavearns

    ltotpen  <- 0.0
    ppregrid <- rep( 0, ngpm ) #pension grid before tax
    pgrid    <- ppregrid #after tax

    for ( im in 1:ngpm ){
      lmearns <- mgrid[Twork,im]

      if (lmearns<= ssbend1) {
          ppregrid[im] = 0.9*lmearns
      }else if (lmearns<= ssbend2) {
          ppregrid[im] = 0.9*ssbend1 + 0.32*(lmearns - ssbend1)
      }else{
          ppregrid[im] = 0.9*ssbend1 + 0.32*(ssbend2-ssbend1) + 0.15*(lmearns - ssbend2)
      }

      ppregrid[im] = ppregrid[im]*lsspar
      pgrid[im] = ppregrid[im] - FnTax( ppregrid[im]*0.85, p )
    }

    #match benefits based on guy with average AIME
    lmearns =  sum( pmin(avearnspre,pencap) )/Twork
    if (lmearns<= ssbend1) {
        lmpen = 0.9*lmearns
    }else if (lmearns<= ssbend2) {
        lmpen = 0.9*ssbend1 + 0.32*(lmearns - ssbend1)
    }else {
        lmpen = 0.9*ssbend1 + 0.32*(ssbend2-ssbend1) + 0.15*(lmearns - ssbend2)
    }
    lmpen = lmpen*lsspar

    if (Display==1) cat('SS Ben for Mean AIME/ Av Pre-tax Earns: ', (lmpen/(sum(avearnspre)/Twork))*100, '%\n')

    if(moment){
      rr = lmpen/(sum(avearnspre)/Twork) - targetSSAvReplacement
    }else{
      rr = list(ppregrid = ppregrid,
                pgrid    = pgrid   )
    }

  })
  return(res)
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

#####################################################
comp.ngpz <- function(iz, p, model, muc, it){
  res <- with( c(p,model), {
    #source('inc.modelsolver_abbg_mpi.r')
    #require(Hmisc)

    emuc <- rep(0,ngpa)
    lcon  <- array( 0, dim=c(ngpa,ngpm,ngpe) )
    lass  <- lcon
    lmuc  <- lcon
    lcon1 <- lcon
    lass1 <- lcon
    lmuc1 <- lcon

    for( ie in 1:ngpe ){
      for( im in 1:ngpm ){
        #solve on tt+1 grid
        if (it==Twork) {
          lmuc1[,im,ie] = bet*Rnet*mucret[1,,im]      #euler equation
        }else{
          emuc[] = 0.0
          for( iz2 in 1:ngpz ){
            for( ie2 in 1:ngpe ){
              #find two probabilities and average over:
              lnextm = ( it*mgrid[it,im] + min(ypregrid[it+1,iz2,ie2],pencap) )/(it+1)
              im2 = FindLinProb1(lnextm, mgrid[it+1,])

              emuc = emuc + im2[2]*muc[it+1,,im2[1]  ,iz2,ie2] *ztrans[it,iz,iz2]*edist +
                            im2[3]*muc[it+1,,im2[1]+1,iz2,ie2] *ztrans[it,iz,iz2]*edist

              if( any(!is.finite(emuc)) ) cat('iz2= ',iz2,'ie2= ',ie2,'\n')
            }
          }
          lmuc1[,im,ie] =  bet*Rnet*emuc                #euler equation
        }

        lcon1[,im,ie] = lmuc1[,im,ie]^(-1.0/gam)

        if ( any(!is.finite(lcon1[,im,ie])) ) cat( 'nan encountered\n')

        if (it==Twork) {
          lass1[,im,ie] = (lcon1[,im,ie] +agridret[1,,im] - ygrid[it,iz,ie] - tran[it,,iz,ie])/Rnet
        }else{
          lass1[,im,ie] = (lcon1[,im,ie] +agrid[it+1,] - ygrid[it,iz,ie] - tran[it,,iz,ie])/Rnet
        }

        #deal with borrowing limits
        if ( min(agrid[it,]) >= lass1[1,im,ie] ) {        #BL does not bind anywhere
          BLbind = 0
        }else{
          #find point in tt grid where BL starts to bind
          BLbind = which.max( agrid[it,][ agrid[it,] < lass1[1,im,ie] ] )
          if (it==Twork) {
            lass[1:BLbind,im,ie] = agridret[1,1,im]
          }else{
            lass[1:BLbind,im,ie] = agrid[it+1,1]
          }
          lcon[1:BLbind,im,ie] = xgrid[it,1:BLbind,iz,ie]- lass[1:BLbind,im,ie]
          lmuc[1:BLbind,im,ie] = lcon[1:BLbind,im,ie]^(-gam)
        }

        #interpolate muc1 as fun of ass1 to get muc where BL does not bind
        lcon[(BLbind+1):ngpa,im,ie] = approxExtrap( lass1[,im,ie], lcon1[,im,ie], agrid[it,(BLbind+1):ngpa] )$y
        lmuc[(BLbind+1):ngpa,im,ie] = lcon[(BLbind+1):ngpa,im,ie]^(-gam)
        lass[(BLbind+1):ngpa,im,ie] = xgrid[it,(BLbind+1):ngpa,iz,ie] - lcon[(BLbind+1):ngpa,im,ie]  #budget constraint
      }
    }

    list(
      lcon      = lcon     ,
      lass      = lass     ,
      lmuc      = lmuc     ,
      lcon1     = lcon1    ,
      lass1     = lass1    ,
      lmuc1     = lmuc1
    )
  })
  return(res)
}
