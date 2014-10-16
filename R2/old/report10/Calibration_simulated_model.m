clear all
clc;

close all

% Canonical earnings process 
load('rw.mat')

N=50000;
T=35;

delta=.2;
tau=.15;
sig=sqrt(.15+.01*((1:1:T)-1));

% check
%std(zsim)./sig

eta_var=zeros(N,T);
eta_var(:,1)=zsim(:,1);
thresh=zeros(T,1);
sig_v=zeros(T-1,1);


for tt=2:T

    % lower threshold; higher threshold is symmetric 
    thresh(tt-1,1)=quantile(eta_var(:,tt-1),1-tau);

    % Var of mean(eta_t|eta_t-1)

    index=eta_var(:,tt-1).*((eta_var(:,tt-1)>=-thresh(tt-1,1)).*(eta_var(:,tt-1)<=thresh(tt-1,1)))+(1-delta*tau)*eta_var(:,tt-1).*((eta_var(:,tt-1)<-thresh(tt-1,1))+(eta_var(:,tt-1)>thresh(tt-1,1)));
    %var(index)

    % Mean of var(eta_t|eta_t-1) (net of the variance of V_t)

    index2=tau*(1-tau)*delta^2*(eta_var(:,tt-1).^2).*((eta_var(:,tt-1)<-thresh(tt-1,1))+(eta_var(:,tt-1)>thresh(tt-1,1)));
    %mean(index2)

    sig_v(tt-1)=sqrt(sig(tt).^2-var(index)-mean(index2));

    V=sig_v(tt-1)*randn(N,1);
    indic=(eta_var(:,tt-1)<-sig(tt-1)*norminv(1-tau)).*(V>sig_v(tt-1)*norminv(1-tau))+(eta_var(:,tt-1)>sig(tt-1)*norminv(1-tau)).*(V<-sig_v(tt-1)*norminv(1-tau));
    rho=1-delta*indic;
    eta_var(:,tt)=rho.*eta_var(:,tt-1)+V;

end
