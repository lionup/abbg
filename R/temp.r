	pind = array( 1:p$ngpp, dim=c(p$ngpm, p$nbin, p$neps) )
	pgrid = array(0,dim=c(p$Tret,p$ngpp)) 
	ip=1
	for(im in 1:p$ngpm){
		for(iz in 1:p$nbin){
			for(ie in 1:p$neps){
		    pgrid[,ip] = rep(mgrid[p$nage,im]*0.5,p$Tret)
		    ip=ip+1
		  }  
	  }  
	}

	agrid = array(0, p$ngpa)
	agrid[1] = p$borrowlim
	for(ia in 2:p$ngpa){
	  agrid[ia] = agrid[ia-1] + exp(p$pexpgrid*(ia-1))
  }
  ltemp1 = (p$amax-agrid[1])/(agrid[p$ngpa]-agrid[1])
  agrid = agrid[1] + ltemp1*(agrid-agrid[1])

	tran = array(0,dim=c(p$nage, p$ngpa, p$nbin, p$neps) )
	xgrid = tran
	tranret = array(0,dim=c(p$Tret, p$ngpa, p$ngpp) )
	xgridret = tranret
	#Cash in Hand
	for(ia in 1:p$ngpa){
	  for(iz in 1:p$nbin){
	    for(ie in 1:p$neps){
	      for( it in 1:p$nage){
	        tran[it,ia,iz,ie] = max(p$cfloor-(p$R*agrid[ia]+ygrid[it,iz,ie]-agrid[1]),0)
	        xgrid[it,ia,iz,ie]= p$R*agrid[ia] + ygrid[it,iz,ie] + tran[it,ia,iz,ie]        
	      }
	  	}
	  }

    for(ip in 1:p$ngpp){
      for(it in 1:p$Tret){
        tranret[it,ia,ip] = max(p$cfloor-(p$R*agrid[ia]+pgrid[it,ip]-agrid[1]),0)
        xgridret[it,ia,ip]= p$R*agrid[ia] + pgrid[it,ip] +tranret[it,ia,ip]
      }
    }
	}

smalli=18
smallj=1
for (j in 1:p$nsim){
	everyi = which(is.na(csim[j,]))[1]
	everyi = ifelse(is.na(everyi),32,everyi)
	if(asim[j,1]==0 & everyi<smalli){
		smalli=everyi
		smallj=j
	}
}

mgrid = array(0,dim=c(p$nage,p$ngpm))        
	#equally spaced between 5th perc and median, and median and 95th perc
	for (it in 1:p$nage){
    lysim = quantile( yavsim[,it], c(0.05,0.5,0.95), names=F, na.rm = T )  
    midm = (1+p$ngpm)/2
    mgrid[it,1]    = lysim[1] 
    mgrid[it,midm] = lysim[2] 
    mgrid[it,p$ngpm] = lysim[3] 

    lwidth1 = (mgrid[it,midm]  - mgrid[it,1])   /(midm-1)
    lwidth2 = (mgrid[it,p$ngpm]- mgrid[it,midm])/(midm-1)
    
    for(im in 2:(midm-1)){
      mgrid[it,im] = mgrid[it,im-1] +lwidth1
    }
    for( im in (midm+1):(p$ngpm-1) ) {
      mgrid[it,im] = mgrid[it,im-1] +lwidth2
    }
	}

		for ( l in (p$nage-1):1 ){
		for ( im in 1:p$ngpm){	
			for ( e in 1:p$nbin){
				Vap=0
    		for ( i in 1:p$neps ){
    			for ( j in (mineta[l+1,e]+1):(maxeta[l+1,e]-1) ){
    				lnextm = (l*mgrid[l,im] + ygrid[l+1,j,i])/(l+1)
    				im2 = FindLinProb1(lnextm,mgrid[l+1,])
						Vap  = Vap + epsprob[i] * etaprob[l+1,e,j] * 
							( GothVP(p, agrid, ygrid[l+1,j,i], C[l+1,im2[1],  j,], M[l+1,im2[1],  j,]) * im2[2] + 
							  GothVP(p, agrid, ygrid[l+1,j,i], C[l+1,im2[1]+1,j,], M[l+1,im2[1]+1,j,]) * im2[3] )  # Gothic Va prime	
					}
				}			
				ChiVec = (p$R * p$beta * Vap)^(-1/p$rho) # inverse Euler equation
	 			MuVec  = agrid+ChiVec
	  		M[l,im,e,]  = c(0, MuVec)    # Matrix of interpolation data
	  		C[l,im,e,]  = c(0,ChiVec)          # Matrix of interpolation data
			}
		}
	}		

	model= list(
	  ygrid = ygrid,
	  mgrid = mgrid,
	  pgrid = pgrid,
    agrid = agrid,
	  Mret  = Mret,
	  Cret  = Cret,  
    M    = M,
    C    = C) 
