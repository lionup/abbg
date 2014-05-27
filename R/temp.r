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

