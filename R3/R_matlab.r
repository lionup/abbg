rm(list=ls())
require(R.matlab)
data <- readMat('~/git/abbg/R3/StephaneNew/data_hermite_cons2.mat')
save(data, file='~/git/abbg/R3/mat_new.dat')

load('rw_nbl_copynl.dat')
with(moments,writeMat('rw_nbl_copynl.mat',zsim=zsim,esim=esim,asim=asim,csim=csim,ysim=ysim,ypresim=ypresim))

load('nl_nbl_eps80.dat')
with(moments,writeMat('nl_nbl_eps80.mat',zsim=zsim,esim=esim,asim=asim,csim=csim,ysim=ysim,ypresim=ypresim))

load('persis_nl_nbl_eps80_parallel.dat')
writeMat('persis_nl_nbl_eps80_parallel.mat',persis=t(persis))

#get the parameters from nl simulations
load('~/git/abbg/R3/eta.dat')
p$Vz0       =  varzapprox[1]
p$Veta_rho1 =  mean(varzapprox[-1] - varzapprox[-18]) #0.01326759  #0.0183   #Veta if rho==1

load('egrid_nl.dat')
aa <- rep(0,18)
for(it in 1:18){
  aa[it] <- sum(egrid[it,]^2 * 1/50) - sum(egrid[it,] * 1/50)^2
}
p$Veps      =  mean(aa)