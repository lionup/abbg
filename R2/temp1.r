grida<-function(allow, p, maxincome){

  xa <- array( 0, dim=c(p$nT+1, p$na, p$nt+p$ntr) ) 

  amid  <- p$am * maxincome #second highest bound on asset grid
  ahigh <- p$ah * maxincome  #upper bound of asset grid, for extreme values
 
  xa[,1,] <- array( allow,dim=c(p$nT+1, p$nt+p$ntr) ) + p$minass
  xa[,2,] <- xa[,1,] + 1e-8
  xa[,3,] <- xa[,2,] + 1e-6 
  xa[,4,] <- xa[,3,] + 1e-4
  xa[,5,] <- xa[,4,] + 1e-3 
 
  for (j in 1:(p$nt+p$ntr)){ 

    adif <- (amid[j] - allow) / 3
    # points 4 to p$na1
    const <- ( adif - (xa[,5,j] - allow) ) / (p$na1-5)
    for (i in 6:p$na1) {
       xa[,i,j] = xa[,5,j] + const * (i - 5)
    }

    # points p$na1+1 to p$na1+p$na2
    const = adif / p$na2
    for (i in 1:p$na2) {
      xa[,p$na1+i,j] = xa[,p$na1,j] + const * i 
    }
        
    # points p$na1+p$na2+1 to p$na1+p$na2+p$na3
    const = adif / p$na3
    for (i in 1:p$na3) {
      xa[,p$na1+p$na2+i,j] = xa[,p$na1+p$na2,j] + const * i 
    }
    
    # points p$na1+p$na2+p$na3+1 to 100
    const = (ahigh[j] - amid[j]) / p$na4
    for (i in 1:p$na4) {
      xa[,p$na1+p$na2+p$na3+i,j] = xa[,p$na1+p$na2+p$na3,j] + const * i 
    }
  }  
  
  return(xa)
} 