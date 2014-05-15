
clear all
clc;

close all

load('data_hermite.mat');

Resqtrue=Resqfinal;

Resqtrue_e0=Resqfinal_e0;

Resqtrue_eps=Resqfinal_eps;

b1true=b1;
bLtrue=bL;
b1true_e0=b1_e0;
bLtrue_e0=bL_e0;
b1true_eps=b1_eps;
bLtrue_eps=bL_eps;

% Draws from the prior distribution of eta's and epsilon's

% seeds
rng('shuffle')

Mateta_true=zeros(N,T);

% Proposal, eta_1
V_draw=unifrnd(0,1,N,1);

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

% Proposal, eta_t
for tt=1:T-1
    Mat=zeros(N,(K1+1)*(K2+1));
    for kk1=0:K1
        for kk2=0:K2            
            Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true(:,tt)-meanY)/stdY).*hermite(kk2,(AGE(:,tt+1)-meanAGE)/stdAGE);
        end
    end
    V_draw=unifrnd(0,1,N,1);
    %First quantile
    Mateta_true(:,tt+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1));
    for jtau=2:Ntau
        Mateta_true(:,tt+1)=Mateta_true(:,tt+1)+...
            ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
            (Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
            (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    end
    %Last quantile.
    Mateta_true(:,tt+1)=Mateta_true(:,tt+1)+(Mat*Resqtrue(:,Ntau)).*...
        (V_draw>Vectau(Ntau));
    
    % Laplace tails
    Mateta_true(:,tt+1)=Mateta_true(:,tt+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
end



Mateps_true=zeros(N,T);

for tt=1:T
    % Proposal, eta_0
    V_draw=unifrnd(0,1,N,1);
    
    %First quantile
    Mateps_true(:,tt)=(MatAGE_t((tt-1)*N+1:tt*N,:)*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
    for jtau=2:Ntau
        Mateps_true(:,tt)=Mateps_true(:,tt)+((MatAGE_t((tt-1)*N+1:tt*N,:)*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+MatAGE_t((tt-1)*N+1:tt*N,:)*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    end
    %Last quantile.
    Mateps_true(:,tt)=Mateps_true(:,tt)+(MatAGE_t((tt-1)*N+1:tt*N,:)*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
    
    % Laplace tails
    Mateps_true(:,tt)=Mateps_true(:,tt)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
end

Ytilde=Mateta_true+Mateps_true;

% Persistence in the data

Vect=Y(:,1:T-1);
Vect_dep=Y(:,2:T);

Mat1=[];
for kk1=0:K1
    
        Mat1=[Mat1 hermite(kk1,(Vect(:)-meanY)/stdY)];
    
end

ResP_data=zeros((K1+1),Ntau);

for jtau=1:Ntau
    
    tau=Vectau(jtau);
    
    ResP_data(:,jtau)=rq(Mat1,Vect_dep(:),tau);
    
end

Mat2=zeros(N*(T-1),1);
for kk1=1:K1
    
        Mat2=[Mat2 kk1*hermite(kk1-1,(Vect(:)-meanY)/stdY)./stdY];
    
end

mean(Mat2*ResP_data)

Vect=quantile(Vect(:),Vectau);
Mat3=zeros(Ntau,1);
for kk1=1:K1
    
        Mat3=[Mat3 kk1*hermite(kk1-1,(Vect(:)-meanY)/stdY)./stdY];
   
end

Mat3*ResP_data
    
% Persistence in the simulated data

VectS=Ytilde(:,1:T-1);
VectS_dep=Ytilde(:,2:T);

MatS1=[];
for kk1=0:K1
    
        MatS1=[MatS1 hermite(kk1,(VectS(:)-meanY)/stdY)];
    
end

ResP_Sdata=zeros((K1+1),Ntau);

for jtau=1:Ntau
    
    tau=Vectau(jtau);
    
    ResP_Sdata(:,jtau)=rq(MatS1,VectS_dep(:),tau);
    
end

MatS2=zeros(N*(T-1),1);
for kk1=1:K1
    
        MatS2=[MatS2 kk1*hermite(kk1-1,(VectS(:)-meanY)/stdY)./stdY];
    
end

mean(MatS2*ResP_Sdata)

VectS=quantile(VectS(:),Vectau);
MatS3=zeros(Ntau,1);
for kk1=1:K1
    
        MatS3=[MatS3 kk1*hermite(kk1-1,(VectS(:)-meanY)/stdY)./stdY];
   
end

MatS3*ResP_Sdata

figure
surf(Vectau,Vectau,Mat3*ResP_data)
axis([0 1 0 1 0 1.2])
xlabel('percentile \tau_{init}')
ylabel('percentile \tau_{shock}')
zlabel('persistence')
figure
surf(Vectau,Vectau,MatS3*ResP_Sdata)
axis([0 1 0 1 0 1.2])
xlabel('percentile \tau_{init}')
ylabel('percentile \tau_{shock}')
zlabel('persistence')


% Persistence of eta

Vect=Mateta_true(:,1:T-1);
    

Vect=quantile(Vect(:),Vectau);
age_ref=meanAGE;
Mat=zeros(Ntau,K2+1);
for kk1=1:K1
    for kk2=0:K2
        Mat=[Mat kk1*hermite(kk1-1,(Vect(:)-meanY)/stdY)./stdY.*hermite(kk2,(age_ref-meanAGE)/stdAGE)];
    end
end


Mat*Resqtrue
figure
surf(Vectau,Vectau,Mat*Resqtrue)
axis([0 1 0 1 0 1.2])
xlabel('percentile \tau_{shock}')
ylabel('percentile \tau_{init}')
zlabel('persistence')

% Densities
[f_eta,xi]=ksdensity(Mateta_true(:));
figure
plot(xi,f_eta)
axis([-2 2 0 1.4])


[f_eps,xi2]=ksdensity(Mateps_true(:));
figure
plot(xi2,f_eps)
axis([-1 1 0 7])
