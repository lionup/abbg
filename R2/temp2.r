# =====================================================================
grida5 <- function(allow, p, maxincome){

  xa <- array( 0, dim=c(p$na, p$nT+1, p$nw, p$nt+p$ntr) ) 

  amid  <- array( rep(p$am*maxincome, each=(p$nT+1)*p$nw), dim=c(p$nT+1, p$nw, p$nt+p$ntr) ) #second highest bound on asset grid
  ahigh <- array( rep(p$ah*maxincome, each=(p$nT+1)*p$nw), dim=c(p$nT+1, p$nw, p$nt+p$ntr) )  #upper bound of asset grid, for extreme values
 
  xa[1,,,] <- allow + p$minass
  xa[2,,,] <- xa[1,,,] + 1e-8
  xa[3,,,] <- xa[2,,,] + 1e-6 
  xa[4,,,] <- xa[3,,,] + 1e-4
  xa[5,,,] <- xa[4,,,] + 1e-3 
 
  adif <- (amid - allow) / 3
  # points 6 to p$na1
  const <- ( adif - (xa[5,,,] - allow) ) / (p$na1-5)
  for (i in 6:p$na1) {
     xa[i,,,] = xa[5,,,] + const * (i - 5)
  }

  # points p$na1+1 to p$na1+p$na2
  const = adif / p$na2
  for (i in 1:p$na2) {
    xa[p$na1+i,,,] = xa[p$na1,,,] + const * i 
  }
      
  # points p$na1+p$na2+1 to p$na1+p$na2+p$na3
  const = adif / p$na3
  for (i in 1:p$na3) {
    xa[p$na1+p$na2+i,,,] = xa[p$na1+p$na2,,,] + const * i 
  }
  
  # points p$na1+p$na2+p$na3+1 to 100
  const = (ahigh - amid) / p$na4
  for (i in 1:p$na4) {
    xa[p$na1+p$na2+p$na3+i,,,] = xa[p$na1+p$na2+p$na3,,,] + const * i 
  }
    
  return(xa)
} 