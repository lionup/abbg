FnGridTrans <- function(lx, p, moment=TRUE){
  res <- with(p,{ 
    legrid <- rep(0,ngpe)
    ledist <- legrid

    # get boundaries and fill in with equally spaced points
    legrid[1]    = -lx*sqrt(Veps)
    legrid[ngpe] = lx*sqrt(Veps)
    lwidth = (legrid[ngpe]-legrid[1])/(ngpe-1)
    for(ij in 2:(ngpe-1)){
      legrid[ij] = legrid[1] + lwidth*(ij-1)
    }

    # fill in probabilities using normal distribution
    ledist[1] = pnorm( legrid[1]+0.5*lwidth, sd=sqrt(Veps) )  
    for( ij in 2:(ngpe-1) ){
      ltemp1 = pnorm( legrid[ij]+0.5*lwidth, sd=sqrt(Veps) )
      ltemp3 = pnorm( legrid[ij]-0.5*lwidth, sd=sqrt(Veps) )
      ledist[ij] = ltemp1 - ltemp3
    }
    ledist[ngpe] = pnorm( legrid[ngpe]-0.5*lwidth, sd=sqrt(Veps), lower.tail=FALSE )
    ledist = ledist/sum(ledist)

    #find variance
    lvar = sum(legrid^2 * ledist) - sum(legrid*ledist)^2

    if(moment){
      #rr = (Veps -lvar)^2
      rr = Veps -lvar
    }else{
      edist <- array( rep(ledist,each=Twork), dim=c(Twork,ngpe) )
      egrid <- array( rep(legrid,each=Twork), dim=c(Twork,ngpe) )
      rr=list(edist=edist, egrid=egrid)
    }
  })     

  return(res)
}

#####################################################
FnGridPerm <- function(lx, p, varz, Vetavec, moment=TRUE){
  res <- with(p,{ 
    lwidth <- rep(0,Twork)
    lvar <- lwidth
    zgrid    <- array( 0, dim=c(Twork,ngpz) ) #permanent component
    zdist  <- zgrid
    ztrans <- array( 0, dim=c(Twork-1,ngpz,ngpz) )
    
    lc <- sqrt(varz)*qnorm(1-tao) #c
    b <- 0.1 * qnorm(1-tao)
    rho3 <- rep(0,ngpz)

    # get boundaries and fill in with equally spaced points
    for(it in 1:Twork){
      zgrid[it,1]    = -lx*sqrt(varz[it])
      zgrid[it,ngpz] = lx*sqrt(varz[it])
      lwidth[it] = (zgrid[it,ngpz]-zgrid[it,1])/(ngpz-1)
      for( iz1 in 2:(ngpz-1)){
        zgrid[it,iz1] = zgrid[it,1] + lwidth[it]*(iz1-1)
      }        
    }

    # fill in transition matrix using normal distribution
    for( it in 1:(Twork-1) ){
      iab = zgrid[it,] >  lc[it]
      ibe = zgrid[it,] < -lc[it]
      ino = !(iAbove | iBelow)

      #normal grids
      for(iz1 in which(ino)){
        ztrans[it,iz1,1] <- pnorm( zgrid[it+1,1]+0.5*lwidth[it+1]-rho*zgrid[it,iz1], sd=sqrt(Vetavec[it]) ) 
        for( iz2 in 2:(ngpz-1) ){
          ltemp1 = pnorm( zgrid[it+1,iz2]+0.5*lwidth[it+1]-rho*zgrid[it,iz1], sd=sqrt(Vetavec[it]) )
          ltemp3 = pnorm( zgrid[it+1,iz2]-0.5*lwidth[it+1]-rho*zgrid[it,iz1], sd=sqrt(Vetavec[it]) )
          ztrans[it,iz1,iz2] = ltemp1 - ltemp3
        }
        ztrans[it,iz1,ngpz] = pnorm( zgrid[it+1,ngpz]-0.5*lwidth[it+1]-rho*zgrid[it,iz1], 
          sd=sqrt(Vetavec[it]), lower.tail=FALSE )
        ztrans[it,iz1,] = ztrans[it,iz1,]/sum(ztrans[it,iz1,])
      }  

      #above c
      for(iz1 in which(iab)){
        #find b
        inl = zgrid[it+1,] < (zgrid[it,iz1]-b)
        rho3[] <- rho 
        rho3[inl] <- rho2

        ztrans[it,iz1,1] <- pnorm( zgrid[it+1,1]+0.5*lwidth[it+1]-rho3[1]*zgrid[it,iz1], sd=sqrt(Vetavec[it]) ) 
        for( iz2 in 2:(ngpz-1) ){
          ltemp1 = pnorm( zgrid[it+1,iz2]+0.5*lwidth[it+1]-rho3[iz2]*zgrid[it,iz1], sd=sqrt(Vetavec[it]) )
          ltemp3 = pnorm( zgrid[it+1,iz2]-0.5*lwidth[it+1]-rho3[iz2]*zgrid[it,iz1], sd=sqrt(Vetavec[it]) )
          ztrans[it,iz1,iz2] = ltemp1 - ltemp3
        }
        ztrans[it,iz1,ngpz] = pnorm( zgrid[it+1,ngpz]-0.5*lwidth[it+1]-rho3[ngpz]*zgrid[it,iz1], 
          sd=sqrt(Vetavec[it]), lower.tail=FALSE )
        ztrans[it,iz1,] = ztrans[it,iz1,]/sum(ztrans[it,iz1,])
      }  

      #below c
      for(iz1 in which(ibe)){
        #find b
        irw = zgrid[it+1,] <= (zgrid[it,iz1]+b) #rw points
        inl = zgrid[it+1,] > (rho2*zgrid[it,iz1]+b) #nl points
        i0  = !(irw|inl) #won't transit into
        ztrans[it,iz1,i0] = 0

        #give autocorre
        rho3[] <- rho 
        rho3[inl] <- rho2

        for(iz2 in which(irw)){
          
        }
        ztrans[it,iz1,1] <- pnorm( zgrid[it+1,1]+0.5*lwidth[it+1]-rho3[1]*zgrid[it,iz1], sd=sqrt(Vetavec[it]) ) 
        for( iz2 in 2:(ngpz-1) ){
          ltemp1 = pnorm( zgrid[it+1,iz2]+0.5*lwidth[it+1]-rho3[iz2]*zgrid[it,iz1], sd=sqrt(Vetavec[it]) )
          ltemp3 = pnorm( zgrid[it+1,iz2]-0.5*lwidth[it+1]-rho3[iz2]*zgrid[it,iz1], sd=sqrt(Vetavec[it]) )
          ztrans[it,iz1,iz2] = ltemp1 - ltemp3
        }
        ztrans[it,iz1,ngpz] = pnorm( zgrid[it+1,ngpz]-0.5*lwidth[it+1]-rho3[ngpz]*zgrid[it,iz1], 
          sd=sqrt(Vetavec[it]), lower.tail=FALSE )
        ztrans[it,iz1,] = ztrans[it,iz1,]/sum(ztrans[it,iz1,])
      }  
    }

    #find distribution at first period
    zdist[1,1] = pnorm( zgrid[1,1]+0.5*lwidth[1], sd=sqrt(Vz0) )  
    for (iz1 in 2:(ngpz-1) ){
      ltemp1 = pnorm( zgrid[1,iz1]+0.5*lwidth[1], sd=sqrt(Vz0) )
      ltemp3 = pnorm( zgrid[1,iz1]-0.5*lwidth[1], sd=sqrt(Vz0) )
      zdist[1,iz1] = ltemp1 - ltemp3
    }
    zdist[1,ngpz] = pnorm( zgrid[1,ngpz]-0.5*lwidth[1], sd=sqrt(Vz0), lower.tail=FALSE )
    zdist[1,] = zdist[1,]/sum(zdist[1,])

    #find unconditional distributions
    for (it in 2:Twork) {
      zdist[it,] = zdist[it-1,] %*% ztrans[it-1,,]
      zdist[it,] = zdist[it,] /sum(zdist[it,])
    }

    #find variance
    for (it in 1:Twork) lvar[it] = sum(zgrid[it,]^2 * zdist[it,]) - sum(zgrid[it,] * zdist[it,])^2

    if(moment){
      rr = sum((varz -lvar)^2)
    }else{
      rr=list(zdist=zdist, zgrid=zgrid, ztrans=ztrans, varzapprox=lvar)
    }
  })     
  return(res)
}
#####################################################
FnGrossInc <- function(lx,lnet,p,stax){
  #lx is gross, lnet is net
  lf = (1-p$pentax - p$btax)*lx + p$btax*( lx^(-p$ptax) + stax )^(-1/p$ptax) - lnet
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

    ltottax = 0.0
    ltotlabincpre = 0.0
    #ltotlabincpost = 0.0

    ygrid    <- array( 0, dim=c(Twork,ngpz,ngpe) ) #earnings
    ypregrid <- ygrid 

    for(it in 1:Twork){
      for(iz in 1:ngpz){
        for(ie in 1:ngpe){
          ygrid[it,iz,ie] = exp( kappa[it] + zgrid[it,iz] + egrid[it,ie] )

          #get implied gross income at this point, to use in constructtion of soc sec system
          lygross = ygrid[it,iz,ie]/(1.0-pentax-btax)
          ypregrid[it,iz,ie] <- uniroot(FnGrossInc, c(0, lygross*3), ygrid[it,iz,ie], p, lstax,
            extendInt="yes", tol=1, maxiter=30)$root

          ltotlabincpre = ltotlabincpre + ypregrid[it,iz,ie]*zdist[it,iz]*edist[it,ie]*popsize[it]
          #ltotlabincpost = ltotlabincpost + ygrid[it,iz,ie]*zdist[it,iz]*edist[it,ie]*popsize[it]
          ltottax = ltottax + btax*( ypregrid[it,iz,ie] - (ypregrid[it,iz,ie]^(-ptax) + lstax)^(-1/ptax) ) *
            zdist[it,iz]*edist[it,ie]*popsize[it]

          avearnspre[it] = avearnspre[it] + ypregrid[it,iz,ie]*zdist[it,iz]*edist[it,ie]
          #avearnspre2[it] = avearnspre2[it] + (ypregrid[it,iz,ie]^2)*zdist[it,iz]*edist[it,ie]
          #avearnspost[it] = avearnspost[it] + ygrid[it,iz,ie]*zdist[it,iz]*edist[it,ie]
          #avearnspost2[it] = avearnspost2[it] + (ygrid[it,iz,ie]^2)*zdist[it,iz]*edist[it,ie]
          #avlearnspre[it] = avlearnspre[it] + log(ypregrid[it,iz,ie])*zdist[it,iz]*edist[it,ie]
          #avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it,iz,ie]^2)*zdist[it,iz]*edist[it,ie]
          #avlearnspost[it] = avlearnspost[it] + log(ygrid[it,iz,ie])*zdist[it,iz]*edist[it,ie]
          #avlearnspost2[it] = avlearnspost2[it] + log(ygrid[it,iz,ie]^2)*zdist[it,iz]*edist[it,ie]
        }
      }
    } 

    #varearnspre    = avearnspre2     - avearnspre^2
    #varearnspost   = avearnspost2   - avearnspost^2
    #varlearnspre   = avlearnspre2   - avlearnspre^2
    #varlearnspost  = avlearnspost2 - avlearnspost^2

    if (Display==1) cat('Tax revenue / Pre-tax labor income: ', (ltottax/ltotlabincpre)*100, '%\n')
    
    if(moment){
      rr = ltottax/ltotlabincpre - targetTaxToLabinc
    }else{
      rr = list(ygrid = ygrid,
                ypregrid = ypregrid, 
                avearnspre   = avearnspre)
    }

  })     
  return(res)
}

#####################################################
FnTax <- function(ly,p, stax){
  lf = p$btax*( ly - (ly^(-p$ptax) + stax)^(-1/p$ptax) ) + p$pentax*ly
}

#####################################################
FnSSParam <- function(lsspar, p, avearnspre, mgrid, pencap, stax, moment=TRUE){
  res <- with(p,{
    if (Display==1) cat('SS benefit parameter: ', lsspar, '\n')

    lavearns = sum(avearnspre)/Twork 
    ssbend1 = 0.18*lavearns
    ssbend2 = 1.1*lavearns
     
    ltotpen = 0.0
    ppregrid <- array( 0, dim=c(Tret,ngpm) ) #pension grid
    pgrid    <- ppregrid

    for ( im in 1:ngpm ){  
      for ( it in 1:Tret ){
        lmearns = mgrid[Twork,im]
        
        if (lmearns<= ssbend1) {
            ppregrid[it,im] = 0.9*lmearns
        }else if (lmearns<= ssbend2) {
            ppregrid[it,im] = 0.9*ssbend1 + 0.32*(lmearns - ssbend1)
        }else{
            ppregrid[it,im] = 0.9*ssbend1 + 0.32*(ssbend2-ssbend1) + 0.15*(lmearns - ssbend2)
        }

        ppregrid[it,im] = ppregrid[it,im]*lsspar
        pgrid[it,im] = ppregrid[it,im] - FnTax( ppregrid[it,im]*0.85, p, stax )
      }        
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

