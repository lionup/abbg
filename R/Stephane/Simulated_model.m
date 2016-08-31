clear all
clc;

close all

% Grid of tau's
Ntau=9;
Vectau=(1/(Ntau+1):1/(Ntau+1):Ntau/(Ntau+1))';




% Nonlinear earnings process 
load('nl_nbl.mat')
age_ref=37-24;
Vect=quantile(zsim(:,age_ref),Vectau);

% Average consumption
meanC_nl=zeros(Ntau,1);
meanC_nl(1)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(1)))/sum((zsim(:,age_ref)<=Vect(1)));
for jtau=2:Ntau-1
     meanC_nl(jtau)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)))/sum((zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)));
end
meanC_nl(Ntau)=sum(csim(:,age_ref).*(zsim(:,age_ref)>Vect(Ntau)))/sum((zsim(:,age_ref)>Vect(Ntau)));

% Average consumption by assets
ql_a=quantile(asim(:,age_ref),.20);
qh_a=quantile(asim(:,age_ref),.80);
meanC_a_nl=zeros(Ntau,2);
meanC_a_nl(1,1)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)<=ql_a))/sum((zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)<=ql_a));
for jtau=2:Ntau-1
     meanC_a_nl(jtau,1)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(jtau+1)).*(asim(:,age_ref)<=ql_a).*(zsim(:,age_ref)>Vect(jtau-1)))/sum((zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)).*(asim(:,age_ref)<=ql_a));
end
meanC_a_nl(Ntau,1)=sum(csim(:,age_ref).*(zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)<=ql_a))/sum((zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)<=ql_a));
meanC_a_nl(1,2)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)>qh_a))/sum((zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)>qh_a));
for jtau=2:Ntau-1
     meanC_a_nl(jtau,2)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(jtau+1)).*(asim(:,age_ref)>qh_a).*(zsim(:,age_ref)>Vect(jtau-1)))/sum((zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)).*(asim(:,age_ref)>qh_a));
end
meanC_a_nl(Ntau,2)=sum(csim(:,age_ref).*(zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)>qh_a))/sum((zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)>qh_a));


% Canonical earnings process 
load('rw_nbl.mat')
Vect=quantile(zsim(:,age_ref),Vectau);

% Average consumption
meanC_rw=zeros(Ntau,1);
meanC_rw(1)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(1)))/sum((zsim(:,age_ref)<=Vect(1)));
for jtau=2:Ntau-1
     meanC_rw(jtau)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)))/sum((zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)));
end
meanC_rw(Ntau)=sum(csim(:,age_ref).*(zsim(:,age_ref)>Vect(Ntau)))/sum((zsim(:,age_ref)>Vect(Ntau)));

% Average consumption by assets
ql_a=quantile(asim(:,age_ref),.20);
qh_a=quantile(asim(:,age_ref),.80);
meanC_a_rw=zeros(Ntau,2);
meanC_a_rw(1,1)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)<=ql_a))/sum((zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)<=ql_a));
for jtau=2:Ntau-1
     meanC_a_rw(jtau,1)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(jtau+1)).*(asim(:,age_ref)<=ql_a).*(zsim(:,age_ref)>Vect(jtau-1)))/sum((zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)).*(asim(:,age_ref)<=ql_a));
end
meanC_a_rw(Ntau,1)=sum(csim(:,age_ref).*(zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)<=ql_a))/sum((zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)<=ql_a));
meanC_a_rw(1,2)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)>qh_a))/sum((zsim(:,age_ref)<=Vect(1)).*(asim(:,age_ref)>qh_a));
for jtau=2:Ntau-1
     meanC_a_rw(jtau,2)=sum(csim(:,age_ref).*(zsim(:,age_ref)<=Vect(jtau+1)).*(asim(:,age_ref)>qh_a).*(zsim(:,age_ref)>Vect(jtau-1)))/sum((zsim(:,age_ref)<=Vect(jtau+1)).*(zsim(:,age_ref)>Vect(jtau-1)).*(asim(:,age_ref)>qh_a));
end
meanC_a_rw(Ntau,2)=sum(csim(:,age_ref).*(zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)>qh_a))/sum((zsim(:,age_ref)>Vect(Ntau)).*(asim(:,age_ref)>qh_a));

% Graph average consumption
figure
Dectau=10*Vectau;
plot(Dectau,[meanC_nl],'--','Linewidth',3,'Color','b')
hold on 
plot(Dectau,[meanC_rw],'-','Linewidth',3,'Color','g')
axis([1 9 5000 65000])
xlabel('decile of \eta_{t-1}','FontSize',20)
ylabel('consumption','FontSize',20)
hold off


% % Graph average consumption (log scale)
% figure
% Dectau=10*Vectau;
% semilogy(Dectau,[meanC_nl],'--','Linewidth',3,'Color','b')
% hold on 
% semilogy(Dectau,[meanC_rw],'-','Linewidth',3,'Color','g')
% axis([1 9 9000 65000])
% xlabel('decile of \eta_{t-1}','FontSize',20)
% ylabel('consumption','FontSize',20)
% hold off



% % Graphs average consumption, by assets
% figure
% plot(Vectau,[meanC_a_nl(:,1) meanC_a_rw(:,1)])
% hold on
% plot(Vectau,[meanC_a_nl(:,2) meanC_a_rw(:,2)])

% Life-cycle evolution: mean, variance, quantiles
load('nl_nbl.mat')
vect_meany_nl=mean(ypresim);
vect_vary_nl=var(ypresim);
vect_quanty_nl=quantile(ypresim,(.1:.1:.9));
vect_meaneta_nl=mean(zsim);
vect_vareta_nl=var(zsim);
vect_meanc_nl=mean(csim);
vect_varc_nl=var(csim);
vect_quantc_nl=quantile(csim,(.1:.1:.9));
vect_meana_nl=mean(asim);
vect_vara_nl=var(asim);
N=50000;
vect_ginia_nl=(N+1)/N-2*sum(sort(asim).*(N+1-(1:1:N)'*ones(1,71)))./(N*sum(asim));
vect_quanta_nl=quantile(asim,(.1:.1:.9));
Mat_nl=cov(ypresim(:,1:20));
load('rw_nbl.mat')
vect_meany_rw=mean(ypresim);
vect_vary_rw=var(ypresim);
vect_quanty_rw=quantile(ypresim,(.1:.1:.9));
vect_meaneta_rw=mean(zsim);
vect_vareta_rw=var(zsim);
vect_meanc_rw=mean(csim);
vect_varc_rw=var(csim);
vect_quantc_rw=quantile(csim,(.1:.1:.9));
vect_meana_rw=mean(asim);
vect_vara_rw=var(asim);
vect_ginia_rw=(N+1)/N-2*sum(sort(asim).*(N+1-(1:1:N)'*ones(1,71)))./(N*sum(asim));
vect_quanta_rw=quantile(asim,(.1:.1:.9));
Mat_rw=cov(ypresim(:,1:20));

% figure
% plot((1:1:35),[vect_meany_nl(1:35);vect_meany_rw(1:35)])
% figure
% plot((1:1:35),[vect_vary_nl(1:35);vect_vary_rw(1:35)])
% figure
% plot((1:1:35),[vect_meaneta_nl(1:35);vect_meaneta_rw(1:35)])
% figure
% plot((1:1:35),[vect_vareta_nl(1:35);vect_vareta_rw(1:35)])
% 
% 
% plot(Dectau,[meanC_nl],'--','Linewidth',3,'Color','b')
% hold on 
% plot(Dectau,[meanC_rw],'-','Linewidth',3,'Color','g')
% axis([1 9 5000 65000])
% xlabel('decile of \eta_{t-1}','FontSize',20)
% ylabel('consumption','FontSize',20)
% hold off


figure
plot((25:1:94),vect_meanc_nl(1:70),'--','Linewidth',3,'Color','b');
hold on
plot((25:1:94),vect_meanc_rw(1:70),'-','Linewidth',3,'Color','g');
axis([25 94 20000 34000])
xlabel('age','FontSize',20)
ylabel('consumption','FontSize',20)
hold off



figure
plot((25:1:94),vect_varc_nl(1:70),'--','Linewidth',3,'Color','b')
hold on
plot((25:1:94),vect_varc_rw(1:70),'-','Linewidth',3,'Color','g')
axis([25 94 0 400000000])
xlabel('age','FontSize',20)
ylabel('consumption variance','FontSize',20)
hold off


% figure
% plot((1:1:70),[vect_meana_nl(1:70);vect_meana_rw(1:70)])
% 
% assets variance
figure
plot((25:1:94),vect_vara_nl(1:70),'--','Linewidth',3,'Color','b')
hold on
plot((25:1:94),vect_vara_rw(1:70),'-','Linewidth',3,'Color','g')
axis([25 94 0 45000000000])
xlabel('age','FontSize',20)
ylabel('assets variance','FontSize',20)
hold off

% assets Gini
% figure
% plot((25:1:94),vect_ginia_nl(1:70),'--','Linewidth',3,'Color','b')
% hold on
% plot((25:1:94),vect_ginia_rw(1:70),'-','Linewidth',3,'Color','g')
% axis([25 94 0 45000000000])
% xlabel('age','FontSize',20)
% ylabel('assets Gini','FontSize',20)
% hold off
% 




% figure
% plot((1:1:35),[vect_quanty_nl(:,1:35);vect_quanty_rw(:,1:35)])
% figure
% plot((1:1:35),[vect_quantc_nl(:,1:35);vect_quantc_rw(:,1:35)])
% figure
% plot((1:1:35),[vect_quanta_nl(:,1:35);vect_quanta_rw(:,1:35)])
% figure
% plot((1:1:10),[Mat_nl(1,1:10);Mat_rw(1,1:10)])
% 



% Impulse responses
% age_ref=37-24;
% load('nl_nbl.mat')
% tau_init=.1;
% 
% qq=quantile(zsim(:,age_ref-1),tau_init);
% indic1=zsim(:,age_ref-1)<=qq;
% yy1=ypresim(indic1,:);
% cc1=csim(indic1,:);
% 
% 
% tau_shock=.1;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)<=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_11_nl=mean(yy);
% meanc_11_nl=mean(cc);
% 
% 
% 
% tau_shock=.5;
% qq=quantile(zsim(indic1,age_ref),tau_shock+.05);
% qqb=quantile(zsim(indic1,age_ref),tau_shock-.05);
% indic2=(zsim(indic1,age_ref)<=qq).*(zsim(indic1,age_ref)>=qqb);
% indic2=(indic2>.5);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_15_nl=mean(yy);
% meanc_15_nl=mean(cc);
% 
% 
% tau_shock=.9;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)>=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_19_nl=mean(yy);
% meanc_19_nl=mean(cc);
% 
% 

% 
% 
% figure
% plot((25:1:59),[meany_11_nl(1:35)-meany_15_nl(1:35);meanc_11_nl(1:35)-meanc_15_nl(1:35)])
% figure
% plot((25:1:59),[meany_19_nl(1:35)-meany_15_nl(1:35);meanc_19_nl(1:35)-meanc_15_nl(1:35)])
% 
% 
% age_ref=37-24;
% load('rw_nbl.mat')
% tau_init=.1;
% 
% qq=quantile(zsim(:,age_ref-1),tau_init);
% indic1=zsim(:,age_ref-1)<=qq;
% yy1=ypresim(indic1,:);
% cc1=csim(indic1,:);
% 
% 
% tau_shock=.1;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)<=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_11_rw=mean(yy);
% meanc_11_rw=mean(cc);
% 
% tau_shock=.5;
% qq=quantile(zsim(indic1,age_ref),tau_shock+.05);
% qqb=quantile(zsim(indic1,age_ref),tau_shock-.05);
% indic2=(zsim(indic1,age_ref)<=qq).*(zsim(indic1,age_ref)>=qqb);
% indic2=(indic2>.5);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_15_rw=mean(yy);
% meanc_15_rw=mean(cc);
% 
% tau_shock=.9;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)>=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_19_rw=mean(yy);
% meanc_19_rw=mean(cc);
% 
% figure
% plot((36:1:59),[meany_11_rw(12:35)-meany_15_rw(12:35);meanc_11_rw(12:35)-meanc_15_rw(12:35)])
% figure
% plot((36:1:59),[meany_19_rw(12:35)-meany_15_rw(12:35);meanc_19_rw(12:35)-meanc_15_rw(12:35)])
% 
% % Impulse responses
% age_ref=37-24;
% load('nl_nbl.mat')
% tau_init=.5;
% 
% qq=quantile(zsim(:,age_ref-1),tau_init);
% indic1=zsim(:,age_ref-1)<=qq;
% yy1=ypresim(indic1,:);
% cc1=csim(indic1,:);
% 
% 
% tau_shock=.1;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)<=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_11_nl=mean(yy);
% meanc_11_nl=mean(cc);
% 
% tau_shock=.5;
% qq=quantile(zsim(indic1,age_ref),tau_shock+.05);
% qqb=quantile(zsim(indic1,age_ref),tau_shock-.05);
% indic2=(zsim(indic1,age_ref)<=qq).*(zsim(indic1,age_ref)>=qqb);
% indic2=(indic2>.5);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_15_nl=mean(yy);
% meanc_15_nl=mean(cc);
% 
% tau_shock=.9;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)>=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_19_nl=mean(yy);
% meanc_19_nl=mean(cc);
% 
% 
% figure
% plot((36:1:59),[meany_11_nl(12:35)-meany_15_nl(12:35);meanc_11_nl(12:35)-meanc_15_nl(12:35)])
% figure
% plot((36:1:59),[meany_19_nl(12:35)-meany_15_nl(12:35);meanc_19_nl(12:35)-meanc_15_nl(12:35)])
% 
% 
% age_ref=37-24;
% load('rw_nbl.mat')
% tau_init=.5;
% 
% qq=quantile(zsim(:,age_ref-1),tau_init);
% indic1=zsim(:,age_ref-1)<=qq;
% yy1=ypresim(indic1,:);
% cc1=csim(indic1,:);
% 
% 
% tau_shock=.1;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)<=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_11_rw=mean(yy);
% meanc_11_rw=mean(cc);
% 
% tau_shock=.5;
% qq=quantile(zsim(indic1,age_ref),tau_shock+.05);
% qqb=quantile(zsim(indic1,age_ref),tau_shock-.05);
% indic2=(zsim(indic1,age_ref)<=qq).*(zsim(indic1,age_ref)>=qqb);
% indic2=(indic2>.5);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_15_rw=mean(yy);
% meanc_15_rw=mean(cc);
% 
% tau_shock=.9;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)>=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_19_rw=mean(yy);
% meanc_19_rw=mean(cc);
% 
% figure
% plot((36:1:59),[meany_11_rw(12:35)-meany_15_rw(12:35);meanc_11_rw(12:35)-meanc_15_rw(12:35)])
% figure
% plot((36:1:59),[meany_19_rw(12:35)-meany_15_rw(12:35);meanc_19_rw(12:35)-meanc_15_rw(12:35)])
% 
% 
% % Impulse responses
% age_ref=37-24;
% load('nl_nbl.mat')
% tau_init=.9;
% 
% qq=quantile(zsim(:,age_ref-1),tau_init);
% indic1=zsim(:,age_ref-1)<=qq;
% yy1=ypresim(indic1,:);
% cc1=csim(indic1,:);
% 
% 
% tau_shock=.1;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)<=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_11_nl=mean(yy);
% meanc_11_nl=mean(cc);
% 
% tau_shock=.5;
% qq=quantile(zsim(indic1,age_ref),tau_shock+.05);
% qqb=quantile(zsim(indic1,age_ref),tau_shock-.05);
% indic2=(zsim(indic1,age_ref)<=qq).*(zsim(indic1,age_ref)>=qqb);
% indic2=(indic2>.5);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_15_nl=mean(yy);
% meanc_15_nl=mean(cc);
% 
% tau_shock=.9;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)>=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_19_nl=mean(yy);
% meanc_19_nl=mean(cc);
% 
% figure
% plot((36:1:59),[meany_11_nl(12:35)-meany_15_nl(12:35);meanc_11_nl(12:35)-meanc_15_nl(12:35)])
% figure
% plot((36:1:59),[meany_19_nl(12:35)-meany_15_nl(12:35);meanc_19_nl(12:35)-meanc_15_nl(12:35)])
% 
% 
% age_ref=37-24;
% load('rw_nbl.mat')
% tau_init=.9;
% 
% qq=quantile(zsim(:,age_ref-1),tau_init);
% indic1=zsim(:,age_ref-1)<=qq;
% yy1=ypresim(indic1,:);
% cc1=csim(indic1,:);
% 
% 
% tau_shock=.1;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)<=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_11_rw=mean(yy);
% meanc_11_rw=mean(cc);
% 
% tau_shock=.5;
% qq=quantile(zsim(indic1,age_ref),tau_shock+.05);
% qqb=quantile(zsim(indic1,age_ref),tau_shock-.05);
% indic2=(zsim(indic1,age_ref)<=qq).*(zsim(indic1,age_ref)>=qqb);
% indic2=(indic2>.5);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_15_rw=mean(yy);
% meanc_15_rw=mean(cc);
% 
% tau_shock=.9;
% qq=quantile(zsim(indic1,age_ref),tau_shock);
% indic2=(zsim(indic1,age_ref)>=qq);
% yy=log(yy1(indic2,:));
% cc=log(cc1(indic2,:));
% meany_19_rw=mean(yy);
% meanc_19_rw=mean(cc);
% 
% figure
% plot((36:1:59),[meany_11_rw(12:35)-meany_15_rw(12:35);meanc_11_rw(12:35)-meanc_15_rw(12:35)])
% figure
% plot((36:1:59),[meany_19_rw(12:35)-meany_15_rw(12:35);meanc_19_rw(12:35)-meanc_15_rw(12:35)])
% 
% 
% 
