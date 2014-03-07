
# integrate x  by log normal
F  <- function(x, Sigma) {
	x * dlnorm(x,-1/2*(Sigma)^2,Sigma)
}	

#CRRA utility function     
u <- function(c,Rho){
	if (Rho==1){
	  u=log(c)
	}else{
	  u = (c^(1-Rho))/(1-Rho)
	}
}

#CRRA marginal utility function      
uP <- function(c,Rho){
	if (Rho==1){
	  uP = 1/c
	}else{
    if(c>0){
      uP= c^(-Rho)
    }else{
      uP= 1e9
    }
  }  
}