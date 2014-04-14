rm(list=ls())
#require(R.matlab)
require(EQL)
require(data.table)
setwd('~/git/abbg/R')
#data <- readMat('data_hermite.mat')
load('mat.dat')
#source('inc.modelsolver.r')
attach(data)

Resqtrue = Resqfinal

Resqtrue_e0 = Resqfinal.e0

Resqtrue_eps = Resqfinal.eps

b1true = b1
bLtrue = bL
b1true_e0 = b1.e0
bLtrue_e0 = bL.e0
b1true_eps = b1.eps
bLtrue_eps = bL.eps

set.seed(123)

# Number of individuals
N=1000000
#N=100
aa_ref=30

nage=15

V_draw <- runif(N)  

MatAGE1 <- rep(0, K3+1)
for (kk3 in 0:K3){    
  MatAGE1[kk3+1]=hermite( (aa_ref-meanAGE)/stdAGE, kk3 )
}

Mateta_true = array(0, dim = c(N,nage))
#Mateps_true = Mateta_true
#Y = Mateta_true

Mateta_true[,1] = (MatAGE1 %*% Resqtrue_e0[,1]) * (V_draw <= Vectau[1])
for (jtau in 2:Ntau){
    Mateta_true[,1] = Mateta_true[,1]+((MatAGE1 %*% (Resqtrue_e0[,jtau]-Resqtrue_e0[,jtau-1]))/(Vectau[jtau]-Vectau[jtau-1]) *
        (V_draw-Vectau[jtau-1]) + MatAGE1 %*% Resqtrue_e0[,jtau-1])*(V_draw>Vectau[jtau-1])*(V_draw<=Vectau[jtau])
}
Mateta_true[,1]=Mateta_true[,1]+(MatAGE1 %*% Resqtrue_e0[,Ntau])*(V_draw>Vectau[Ntau])

Mateta_true[,1]=Mateta_true[,1]+( (1/(b1true_e0)*log(V_draw/Vectau[1]))*(V_draw<=Vectau[1]) - 
    (1/bLtrue_e0*log((1-V_draw)/(1-Vectau[Ntau])))*(V_draw>Vectau[Ntau]))

for (jj in 1:nage){
    
    aa=aa_ref+(jj-1)*2
    
	#for (kk3 in 0:K3){    
	#  MatAGE1[kk3+1]=hermite( (aa-meanAGE)/stdAGE, kk3 )
	#}
    
    #MatAGE_t=[];
    #for kk4=0:K4
    #    MatAGE_t=[MatAGE_t hermite(kk4,(aa-meanAGE)/stdAGE)];
    #end
    
    
    
    #V_draw <- runif(N)  
    
    #Mateps_true(:,jj)=(MatAGE_t*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
    #for jtau=2:Ntau
    #    Mateps_true(:,jj)=Mateps_true(:,jj)+((MatAGE_t*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
    #        (V_draw-Vectau(jtau-1))+MatAGE_t*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    #end
    #%Last quantile.
    #Mateps_true(:,jj)=Mateps_true(:,jj)+(MatAGE_t*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
    
    #% Laplace tails
    #Mateps_true(:,jj)=Mateps_true(:,jj)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
    #    -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
    
    
    # Earnings
    
    # Y(:,jj)=Mateta_true(:,jj)+Mateps_true(:,jj);
    
    # Eta
    
    if (jj <= nage-1){ #this period is jj+1
        
        Mat = array(0, dim=c(N,(K1+1)*(K2+1)))
        for (kk1 in 0:K1){
            for (kk2 in 0:K2) {
                Mat[,kk1*(K2+1)+kk2+1]=hermite( (Mateta_true[,jj]-meanY)/stdY, kk1 ) * hermite( ((aa+2)-meanAGE)/stdAGE, kk2 )
            }
        }
        
        V_draw <- runif(N)  
        
        #First quantile
        
        Mateta_true[,jj+1] = ( Mat %*% Resqtrue[,1] ) * ( V_draw <= Vectau[1] )

        for (jtau in 2:Ntau){
            Mateta_true[,jj+1] = Mateta_true[,jj+1] + 
                ( (Mat %*% (Resqtrue[,jtau]-Resqtrue[,jtau-1])) / 
                (Vectau[jtau]-Vectau[jtau-1]) *
                (V_draw - Vectau[jtau-1]) + Mat %*% Resqtrue[,jtau-1] ) * 
                (V_draw > Vectau[jtau-1]) * (V_draw <= Vectau[jtau])
        }
        Mateta_true[,jj+1]=Mateta_true[,jj+1] + (Mat %*% Resqtrue[,Ntau]) * (V_draw>Vectau[Ntau])

        Mateta_true[,jj+1]=Mateta_true[,jj+1] + ( (1 / b1true * log(V_draw/Vectau[1])) * (V_draw <= Vectau[1]) - 
            (1 / bLtrue * log((1-V_draw)/(1-Vectau[Ntau]))) * (V_draw > Vectau[Ntau]) )

        
    }
}

detach(data)

# Quantiles of eta and epsilon, by age
nbin <- 23
xeta <- array( 0, dim=c(nage, nbin) )  #21 bin
veta <- seq(1/(2*nbin), (2*nbin-1)/(2*nbin), l=(2*nbin-1)) #41 nodes
oddnode <- seq(1,(2*nbin-1),2)  #income node is the even node
age = seq(30, (30+(nage-1)*2), 2)

for (i in 1:nage){
    xeta[i,] <- quantile( Mateta_true[,i], veta, names=F, na.rm = T )[oddnode]
}    


#quantile( Mateta_true[,1], veta, names=F, na.rm = T )
#table(cut_number(Mateta_true[,1],21))
#quantile(Mateta_true,(1/24:1/24:23/24))
#quantile(Mateps_true,(1/24:1/24:23/24))

#for (i in 1:nage){
#    longi <- data.table( pid = 1:N, age=age[i], income = Mateta_true[,i])
#    longi[,q:=cut(income,xeta[i,])]
#    longi[,rank:=ceiling(23*rank(income)/N), age]
#    longi[,qn:=as.numeric(cut_number(income, n = 23))]
#}

long <- data.table( pid = 1:N, age=rep(age, each=N), income = c(Mateta_true))
setkey(long, pid, age)
long[, q:=as.numeric(cut_number(income, n = nbin)), age]
longq <- long
longq[,income:= NULL]
wide <- reshape(longq, idvar='pid', timevar='age', direction='wide')

etaprob <- array( 0,dim=c(nage, nbin, nbin) )
etaprob[1,,] <- 1/nbin

trans.matrix <- function(x1, x2, prob=T)
{
    tt <- table( x1, x2 )
    if(prob) tt <- tt / rowSums(tt)
    tt
}

for (t in 1:(nage-1)){
    x1 <- paste('q',age[t],sep='.')
    x2 <- paste('q',age[t]+2,sep='.')
    etaprob[t+1,,] <- trans.matrix(wide[[x1]],wide[[x2]])
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


etacontot <- etaprob
for (i in 1:nage) {
  etacontot[i,,] <- econdCDF(nbin, etaprob[i,,])
}  
etauntot <- etaprob[1,1,]
for(i in 2:nbin){
    etauntot[i] <- etauntot[i-1]  + etauntot[i]  
}

nsim <- 10000
#simulate permanent income
etalong <- rep(0, nage)
etaind <- array(0, dim=c(nsim,nage ))
randeta <-array(0, dim=c(nsim, nage))
visite  <-array(0, dim=c(nage, nbin))

for (i in 1:nage){
    randeta[,i] <- runif(nsim)  #income
}

#loop for different individuals
for (i in 1:nsim){
 
    #loop over life cycle
    for (t in 1:nage) {     

        # find income draw for individual i at period 1
        if(t==1){
            enode <- which( randeta[i,t] <= etauntot )[1]
            visite[t,enode] = visite[t,enode] + 1
        # find income draw for individual i at period t
        }else{
            enode <- which( randeta[i,t] <= etacontot[t,enode,] )[1]
            visite[t,enode] = visite[t,enode] + 1
        }

        etaind[i,t]<-xeta[t,enode]
        etalong[t] <- etalong[t] + xeta[t,enode]
    }       
}   

etamedian <- apply(etaind,2,median)
