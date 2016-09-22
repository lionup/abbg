FnGridTrans <- function(lx, p, vepsi, moment=TRUE){
  res <- with(p,{
    legrid <- rep(0,ngpe)
    ledist <- legrid

    # get boundaries and fill in with equally spaced points
    legrid[1]    = -lx*sqrt(vepsi)
    legrid[ngpe] = lx*sqrt(vepsi)
    lwidth = (legrid[ngpe]-legrid[1])/(ngpe-1)
    for(ij in 2:(ngpe-1)){
      legrid[ij] = legrid[1] + lwidth*(ij-1)
    }

    # fill in probabilities using normal distribution
    ledist[1] = pnorm( legrid[1]+0.5*lwidth, sd=sqrt(vepsi) )
    for( ij in 2:(ngpe-1) ){
      ltemp1 = pnorm( legrid[ij]+0.5*lwidth, sd=sqrt(vepsi) )
      ltemp3 = pnorm( legrid[ij]-0.5*lwidth, sd=sqrt(vepsi) )
      ledist[ij] = ltemp1 - ltemp3
    }
    ledist[ngpe] = pnorm( legrid[ngpe]-0.5*lwidth, sd=sqrt(vepsi), lower.tail=FALSE )
    ledist = ledist/sum(ledist)

    #find variance
    lvar = sum(legrid^2 * ledist) - sum(legrid*ledist)^2

    if(moment){
      rr = vepsi -lvar
    }else{
      rr=list(edisti=ledist, egridi=legrid)
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

    # get boundaries and fill in with equally spaced points
    for(it in 1:Twork){
      zgrid[it,1]    = - lx*sqrt(varz[it])
      zgrid[it,ngpz] = lx*sqrt(varz[it])
      lwidth[it] = (zgrid[it,ngpz]-zgrid[it,1])/(ngpz-1)
      for( iz1 in 2:(ngpz-1)){
        zgrid[it,iz1] = zgrid[it,1] + lwidth[it]*(iz1-1)
      }
    }

    # fill in transition matrix using normal distribution
    for( it in 1:(Twork-1) ){
      for (iz1 in 1:ngpz){
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
    }

    #find distribution at first period
    zdist[1,1] = pnorm( zgrid[1,1]+0.5*lwidth[1], sd=sqrt(varz[1]) )
    for (iz1 in 2:(ngpz-1) ){
      ltemp1 = pnorm( zgrid[1,iz1]+0.5*lwidth[1], sd=sqrt(varz[1]) )
      ltemp3 = pnorm( zgrid[1,iz1]-0.5*lwidth[1], sd=sqrt(varz[1]) )
      zdist[1,iz1] = ltemp1 - ltemp3
    }
    zdist[1,ngpz] = pnorm( zgrid[1,ngpz]-0.5*lwidth[1], sd=sqrt(varz[1]), lower.tail=FALSE )
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

    load('~/git/abbg/R3/difeta.dat')
    #difeta <- rep(0, Twork)

    for(it in 1:Twork){
      for(iz in 1:ngpz){
        for(ie in 1:ngpe){
          ygrid[it,iz,ie] = exp( kappa[it] + zgrid[it,iz] + egrid[it,ie] ) + difeta[it]


          #get implied gross income at this point, to use in constructtion of soc sec system
          lygross <- ygrid[it,iz,ie]/(1.0-pentax-btax)
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
FnTax <- function(ly,p){
  lf = p$btax*( ly - (ly^(-p$ptax) + p$stax)^(-1/p$ptax) ) + p$pentax*ly
}

#####################################################
FnSSParam <- function(lsspar, p, avearnspre, mgrid, moment=TRUE){
  res <- with(p,{
    if (Display==1) cat('SS benefit parameter: ', lsspar, '\n')

    lavearns <- sum(avearnspre)/Twork
    ssbend1  <- 0.18*lavearns
    ssbend2  <- 1.1*lavearns

    ltotpen  <- 0.0
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
        pgrid[it,im] = ppregrid[it,im] - FnTax( ppregrid[it,im]*0.85, p )
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

              emuc = emuc + im2[2]*muc[it+1,,im2[1]  ,iz2,ie2] *ztrans[it,iz,iz2]*edist[it+1,ie2] +
                            im2[3]*muc[it+1,,im2[1]+1,iz2,ie2] *ztrans[it,iz,iz2]*edist[it+1,ie2]

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
