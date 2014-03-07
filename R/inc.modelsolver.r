
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


# Inverse of the CRRA marginal 
nP <- function(c,Rho){
  nP = c^(-1/Rho)
}



% GothVP.m
% Gothic V prime function

function EUP = GothVP(a,rho)

% Another M-file function, Cnextp, is also involved in running this program

global Rhat uP nP ThetaVec ThetaVecProb PermVec PermVecProb Beta R

% EUP is used to take the sum of marginal utilities UP weighted by probabilities
EUP = zeros(size(a)); 

for i=1:length(ThetaVec)
  for j=1:length(PermVec)
    PermVal = PermVec(j);
    ThetaVal = ThetaVec(i);
    EUP     = EUP + Beta.*R...
    .*uP(PermVal.*Cnextp(R.*a./PermVal+ones(1,length(EUP)).*ThetaVal),rho).*ThetaVecProb(i)...
  .*PermVecProb(j);       
  end;
end;