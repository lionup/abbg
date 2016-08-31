
clear all
clc;

global Sigma T sig2eta1 sig2v sig2eps

load('data_hermite_cons2.mat');


% Earnings
Sigma=zeros(T,T);
for tt=1:T
    for ttt=1:T
        Sigma(tt,ttt)=mean(Y(:,tt).*Y(:,ttt))-mean(Y(:,tt))*mean(Y(:,ttt));
    end
end

par0=[.2;.2;.2]+randn(3,1);

[par fval]=fminunc(@covariance_cm,par0);

sig2eta1=par(1)^2;
sig2v=par(2)^2;
sig2eps=par(3)^2;

fval
par



% Covariance matrix
Sigma=zeros(2*T,2*T);
for tt=1:T
    for ttt=1:T
        Sigma(tt,ttt)=mean(Y(:,tt).*Y(:,ttt))-mean(Y(:,tt))*mean(Y(:,ttt));
    end
    for ttt=T+1:2*T
        Sigma(tt,ttt)=mean(Y(:,tt).*C(:,ttt-T))-mean(Y(:,tt))*mean(C(:,ttt-T));
    end
end

for tt=T+1:2*T
    for ttt=1:T
        Sigma(tt,ttt)=mean(C(:,tt-T).*Y(:,ttt))-mean(C(:,tt-T))*mean(Y(:,ttt));
    end
    for ttt=T+1:2*T
        Sigma(tt,ttt)=mean(C(:,tt-T).*C(:,ttt-T))-mean(C(:,tt-T))*mean(C(:,ttt-T));
    end
end


par0=[0;0;.2]+randn(3,1);

[par fval]=fminunc(@covariance_cm_consumption_twostep,par0);

phi_eta=par(1);
phi_eps=par(2);
sig2xi=par(3)^2;

fval
par

Sig=zeros(2*T,2*T);

% (Y,Y)
for tt=1:T
    Sig(tt,tt)=sig2eta1+(tt-1)*sig2v+sig2eps;
end

for tt=1:T
    for ttt=tt+1:T
        Sig(tt,ttt)=sig2eta1+(tt-1)*sig2v;
        Sig(ttt,tt)=sig2eta1+(tt-1)*sig2v;
    end
    
end

% (Y,C)
for tt=1:T
    Sig(tt,T+tt)=phi_eta*(sig2eta1+(tt-1)*sig2v)+phi_eps*sig2eps;
end

for tt=1:T
    for ttt=tt+1:T
        Sig(tt,T+ttt)=phi_eta*(sig2eta1+(tt-1)*sig2v);
        Sig(ttt,T+tt)=phi_eta*(sig2eta1+(tt-1)*sig2v);
    end
    
end

% (C,Y)
for tt=1:T
    Sig(T+tt,tt)=phi_eta*(sig2eta1+(tt-1)*sig2v)+phi_eps*sig2eps;
end

for tt=1:T
    for ttt=tt+1:T
        Sig(T+tt,ttt)=phi_eta*(sig2eta1+(tt-1)*sig2v);
        Sig(T+ttt,tt)=phi_eta*(sig2eta1+(tt-1)*sig2v);
    end
    
end

% (C,C)
for tt=1:T
    Sig(T+tt,T+tt)=phi_eta^2*(sig2eta1+(tt-1)*sig2v)+phi_eps^2*sig2eps+sig2xi;
end

for tt=1:T
    for ttt=tt+1:T
        Sig(T+tt,T+ttt)=phi_eta^2*(sig2eta1+(tt-1)*sig2v);
        Sig(T+ttt,T+tt)=phi_eta^2*(sig2eta1+(tt-1)*sig2v);
    end
    
end

Sigma

Sig

% Consumption response to permanent income shocks

Ntau2=30;
Vectau2=(1/(Ntau2+1):1/(Ntau2+1):Ntau2/(Ntau2+1))';

Vect_impact=phi_eta*sqrt(sig2v)./normpdf(norminv(Vectau2));
plot(Vectau2,Vect_impact)

% Simulation

% Expand the sample
Nsim=20;
N=N*Nsim;


Ytilde=zeros(N,T);
Mateta_true=zeros(N,T);

Mateta_true(:,1)=sqrt(sig2eta1)*randn(N,1);

Mateps_true=sqrt(sig2eps)*randn(N,T);



for tt=1:T
   
    Ytilde(:,tt)=Mateta_true(:,tt)+Mateps_true(:,tt);
    
    
    if tt<=T-1
        Mateta_true(:,tt+1)=Mateta_true(:,tt)+sqrt(sig2v)*randn(N,1);
        
    end
end


% Persistence

K1=3;
meanY=mean(Y(:));
stdY=std(Y(:));

% Grid of tau's
Ntau=11;
Vectau=(1/(Ntau+1):1/(Ntau+1):Ntau/(Ntau+1))';


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
surf(Vectau,Vectau,MatS3*ResP_Sdata)
xlabel('percentile \tau_{shock}','FontSize',14)
ylabel('percentile \tau_{init}','FontSize',14)
zlabel('persistence','FontSize',14)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1.2])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1.2))
