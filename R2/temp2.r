		for( it in Twork:1){
		  if (Display==1) cat( 'Solving for decision rules at age ', it, '\n')             

		  for( iz in 1:ngpz ){
    for( ie in 1:ngpe ){
      for( im in 1:ngpm ){        
        #solve on tt+1 grid
        if (it==Twork) {         
          muc1[it,,im,iz,ie] = bet*Rnet*mucret[1,,im]      #euler equation
        }else{
          emuc[] = 0.0
          for( iz2 in 1:ngpz ){
            for( ie2 in 1:ngpe ){
              #find two probabilities and average over:
              lnextm = ( it*mgrid[it,im] + min(ypregrid[it+1,iz2,ie2],pencap) )/(it+1)
              im2 = FindLinProb1(lnextm, mgrid[it+1,])
              
              emuc = emuc + im2[2]*muc[it+1,,im2[1]  ,iz2,ie2] *ztrans[it,iz,iz2]*edist[it+1,ie2] +
                            im2[3]*muc[it+1,,im2[1]+1,iz2,ie2] *ztrans[it,iz,iz2]*edist[it+1,ie2]                         
            
              if( any(is.na(emuc)) ) cat('iz2= ',iz2,'ie2= ',ie2,'\n')
            }
          }
          muc1[it,,im,iz,ie] =  bet*Rnet*emuc                #euler equation
        }
        
        con1[it,,im,iz,ie] = muc1[it,,im,iz,ie]^(-1.0/gam) 
        
        if ( any(!is.finite(con1[it,,im,iz,ie])) ) cat( 'nan encountered\n')
              
        if (it==Twork) { 
          ass1[it,,im,iz,ie] = (con1[it,,im,iz,ie] +agridret[1,,im] - ygrid[it,iz,ie] - tran[it,,iz,ie])/Rnet  
        }else{
          ass1[it,,im,iz,ie] = (con1[it,,im,iz,ie] +agrid[it+1,] - ygrid[it,iz,ie] - tran[it,,iz,ie])/Rnet
        }
        
        #deal with borrowing limits
        if ( min(agrid[it,]) >= ass1[it,1,im,iz,ie] ) {        #BL does not bind anywhere
          BLbind = 0  
        }else{
          #find point in tt grid where BL starts to bind
          BLbind = which.max( agrid[it,][ agrid[it,] < ass1[it,1,im,iz,ie] ] )
          if (it==Twork) { 
            ass[it,1:BLbind,im,iz,ie] = agridret[1,1,im]
          }else{
            ass[it,1:BLbind,im,iz,ie] = agrid[it+1,1]
          }
          con[it,1:BLbind,im,iz,ie] = xgrid[it,1:BLbind,iz,ie]-ass[it,1:BLbind,im,iz,ie]
          muc[it,1:BLbind,im,iz,ie] = con[it,1:BLbind,im,iz,ie]^(-gam)
        }

        #interpolate muc1 as fun of ass1 to get muc where BL does not bind      
        con[it,(BLbind+1):ngpa,im,iz,ie] = approxExtrap( ass1[it,,im,iz,ie], con1[it,,im,iz,ie], agrid[it,(BLbind+1):ngpa] )$y
        muc[it,(BLbind+1):ngpa,im,iz,ie] = con[it,(BLbind+1):ngpa,im,iz,ie]^(-gam)        
        ass[it,(BLbind+1):ngpa,im,iz,ie] = xgrid[it,(BLbind+1):ngpa,iz,ie] -con[it,(BLbind+1):ngpa,im,iz,ie]  #budget constraint
      }
    }
		  }
		} 