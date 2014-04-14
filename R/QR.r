rm(list=ls())
load(mat.dat)
attach(data)

Resqtrue=Resqfinal
Resqtrue_e0=Resqfinal.e0
Resqtrue_eps=Resqfinal.eps

b1true=b1
bLtrue=bL
b1true_e0=b1.e0
bLtrue_e0=bL.e0
b1true_eps=b1.eps
bLtrue_eps=bL.eps

#Draws from the prior distribution of eta's and epsilon's

Mateta_true=array(0, dim=c(N,T))

# Proposal, eta_1
V_draw=runif(N)

#First quantile
Mateta_true[,1]=(MatAGE1*Resqtrue_e0[,1])*(V_draw<=Vectau(1))
for (jtau in 2:Ntau){
    Mateta_true[,1]=Mateta_true[,1]+((MatAGE1*(Resqtrue_e0(:,jtau)-Resqtrue_e0(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
        (V_draw-Vectau(jtau-1))+MatAGE1*Resqtrue_e0(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau))
}
#Last quantile.
Mateta_true(:,1)=Mateta_true(:,1)+(MatAGE1*Resqtrue_e0(:,Ntau)).*(V_draw>Vectau(Ntau))

# Laplace tails
Mateta_true(:,1)=Mateta_true(:,1)+((1/(b1true_e0)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
    -(1/bLtrue_e0*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)))

# Proposal, eta_t
for tt=1:T-1
    Mat=zeros(N,(K1+1)*(K2+1))
    for kk1=0:K1
        for kk2=0:K2            
            Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true(:,tt)-meanY)/stdY).*hermite(kk2,(AGE(:,tt+1)-meanAGE)/stdAGE)
        end
    end
    V_draw=unifrnd(0,1,N,1)
    %First quantile
    Mateta_true(:,tt+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1))
    for jtau=2:Ntau
        Mateta_true(:,tt+1)=Mateta_true(:,tt+1)+...
            ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
            (Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
            (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau))
    end
    %Last quantile.
    Mateta_true(:,tt+1)=Mateta_true(:,tt+1)+(Mat*Resqtrue(:,Ntau)).*...
        (V_draw>Vectau(Ntau))
    
    % Laplace tails
    Mateta_true(:,tt+1)=Mateta_true(:,tt+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)))
    
end

