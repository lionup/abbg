
clear all
clc;

close all

%%3d plot consumption response conditional on asset and age

load('~/git/abbg/R3/persis_nl_nbl_e50m50_parallel.mat')
load('~/git/abbg/R3/StephaneNew/data_hermite_cons2.mat')

figure
surf(Vectau,Vectau,persis)
xlabel('percentile \tau_{age}','FontSize',14)
ylabel('percentile \tau_{assets}','FontSize',14)
zlabel('consumption response','FontSize',14)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 .8])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:.2:.8))
colormap(jet)


%% 3d plot nonlinear Persistence in simulation
clear
load('~/git/abbg/R/figure/report14/persisinc_cohort30_parallel.mat')
load('~/git/abbg/R/data_hermite.mat')

figure
surf(Vectau,Vectau,persisinc)
xlabel('percentile \tau_{shock}','FontSize',14)
ylabel('percentile \tau_{init}','FontSize',14)
zlabel('persistence','FontSize',14)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1.2])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1.2))
colormap(jet)

%% Persistence in the data
clear
load('~/git/abbg/R/figure/report14/persisinc_data_cohort30_parallel.mat')
load('~/git/abbg/R/data_hermite.mat')

figure
surf(Vectau,Vectau,persisinc)
xlabel('percentile \tau_{shock}','FontSize',14)
ylabel('percentile \tau_{init}','FontSize',14)
zlabel('persistence','FontSize',14)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1.2])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1.2))
colormap(jet)