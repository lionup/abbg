rm(list=ls())
require(R.matlab)
data <- readMat('~/git/abbg/R3/StephaneNew/data_hermite_cons2.mat')
save(data, file='~/git/abbg/R3/mat_new.dat')

load('rw_nbl_full.dat')
with(moments,writeMat('rw_nbl_full.mat',zsim=zsim,esim=esim,asim=asim,csim=csim,ysim=ysim,ypresim=ypresim))

load('nl_nbl_e50m50.dat')
with(moments,writeMat('nl_nbl_e50m50.mat',zsim=zsim,esim=esim,asim=asim,csim=csim,ysim=ysim,ypresim=ypresim))