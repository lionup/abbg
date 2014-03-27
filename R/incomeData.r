rm(list=ls())
#require(R.matlab)
require(EQL)
#require(data.table)
setwd('~/git/abbg/R')
#data <- readMat('data_hermite.mat')
load('mat.dat')
source('inc.modelsolver.r')
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

aa_ref=30

nage=18

V_draw <- runif(N)  

MatAGE1 <- rep(0, K3+1)
for (kk3 in 0:K3){    
  MatAGE1[kk3+1]=hermite( (aa_ref-meanAGE)/stdAGE, kk3 )
}

Mateta_true = array(0, dim = c(N,nage))
Mateps_true = Mateta_true
Y = Mateta_true

Mateta_true[,1] = (MatAGE1 %*% Resqtrue_e0[,1]) * (V_draw <= Vectau[1])
for (jtau in 2:Ntau){
    Mateta_true[,1] = Mateta_true[,1]+((MatAGE1 %*% (Resqtrue_e0[,jtau]-Resqtrue_e0[,jtau-1]))/(Vectau[jtau]-Vectau[jtau-1]) *
        (V_draw-Vectau[jtau-1]) + MatAGE1 %*% Resqtrue_e0[,jtau-1])*(V_draw>Vectau[jtau-1])*(V_draw<=Vectau[jtau])
}
Mateta_true[,1]=Mateta_true[,1]+(MatAGE1 %*% Resqtrue_e0[,Ntau])*(V_draw>Vectau[Ntau])

Mateta_true[,1]=Mateta_true[,1]+( (1/(b1true_e0)*log(V_draw/Vectau[1]))*(V_draw<=Vectau[1]) - 
    (1/bLtrue_e0*log((1-V_draw)/(1-Vectau[Ntau])))*(V_draw>Vectau[Ntau]))

for (jj in 2:nage){
    
    aa=aa_ref+(jj-1)*2
    
		for (kk3 in 0:K3){    
		  MatAGE1[kk3+1]=hermite( (aa-meanAGE)/stdAGE, kk3 )
		}
    
    MatAGE_t=[];
    for kk4=0:K4
        MatAGE_t=[MatAGE_t hermite(kk4,(aa-meanAGE)/stdAGE)];
    end
    
    
    
    V_draw=unifrnd(0,1,N,1);
    
    %First quantile
    Mateps_true(:,jj)=(MatAGE_t*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
    for jtau=2:Ntau
        Mateps_true(:,jj)=Mateps_true(:,jj)+((MatAGE_t*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+MatAGE_t*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    end
    %Last quantile.
    Mateps_true(:,jj)=Mateps_true(:,jj)+(MatAGE_t*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
    
    % Laplace tails
    Mateps_true(:,jj)=Mateps_true(:,jj)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
    
    
    % Earnings
    
    Y(:,jj)=Mateta_true(:,jj)+Mateps_true(:,jj);
    
    
    
    
    % Eta
    
    if jj<=nage-1
        
        Mat=zeros(N,(K1+1)*(K2+1));
        for kk1=0:K1
            for kk2=0:K2
                Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true(:,jj)-meanY)/stdY).*hermite(kk2,((aa+2)-meanAGE)/stdAGE);
            end
        end
        
        V_draw=unifrnd(0,1,N,1);
        
        %First quantile
        
        Mateta_true(:,jj+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1));
        for jtau=2:Ntau
            Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+...
                ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
                (Vectau(jtau)-Vectau(jtau-1)).*...
                (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
                (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
        end
        %Last quantile.
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+(Mat*Resqtrue(:,Ntau)).*...
            (V_draw>Vectau(Ntau));
        
        % Laplace tails
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
            -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
        
    end
end

% Quantiles of eta and epsilon, by age

quantile(Mateta_true,(1/24:1/24:23/24))
quantile(Mateps_true,(1/24:1/24:23/24))
