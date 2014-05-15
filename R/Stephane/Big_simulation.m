
clear all
clc;

close all


load('data_hermite_cons.mat');

Resqtrue=Resqfinal;

Resqtrue_e0=Resqfinal_e0;

Resqtrue_eps=Resqfinal_eps;

b1true=b1;
bLtrue=bL;
b1true_e0=b1_e0;
bLtrue_e0=bL_e0;
b1true_eps=b1_eps;
bLtrue_eps=bL_eps;



% seeds
rng('shuffle')

% Number of individuals
N=1000000;

% starting age, and number of periods (nage*2years)

aa_ref=30;

nage=15;

% First period

% eta_1
V_draw=unifrnd(0,1,N,1);

MatAGE1=[];
for kk3=0:K3
    MatAGE1=[MatAGE1 hermite(kk3,(aa_ref-meanAGE)/stdAGE)];
end

Mateta_true=zeros(N,nage);
Mateps_true=zeros(N,nage);
Y=zeros(N,nage);



%First quantile
Mateta_true(:,1)=(MatAGE1*Resqtrue_e0(:,1)).*(V_draw<=Vectau(1));
for jtau=2:Ntau
    Mateta_true(:,1)=Mateta_true(:,1)+((MatAGE1*(Resqtrue_e0(:,jtau)-Resqtrue_e0(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
        (V_draw-Vectau(jtau-1))+MatAGE1*Resqtrue_e0(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
end
%Last quantile.
Mateta_true(:,1)=Mateta_true(:,1)+(MatAGE1*Resqtrue_e0(:,Ntau)).*(V_draw>Vectau(Ntau));

% Laplace tails
Mateta_true(:,1)=Mateta_true(:,1)+((1/(b1true_e0)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
    -(1/bLtrue_e0*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));

% Subsequent periods


for jj=1:nage
    
    aa=aa_ref+(jj-1)*2;
    
    MatAGE1=[];
    for kk3=0:K3
        MatAGE1=[MatAGE1 hermite(kk3,(aa-meanAGE)/stdAGE)];
    end
    
    % epsilons
    
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
