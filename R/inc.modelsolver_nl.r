    trans.matrix <- function(x1, x2, prob=T)
    {
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