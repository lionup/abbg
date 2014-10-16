#####################################################
comp.ngpz <- function(iz, p, model, muc, it){
  res <- with( c(p,model), {
    #source('inc.modelsolver_abbg_mpi.r')
    require(Hmisc)

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
