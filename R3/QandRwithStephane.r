require(ggplot2)
require(data.table)

load('~/git/abbg/R2/rw_nbl.dat')
moments_switch <- moments

mean_y_switch <- apply(moments_switch$ysim[,1:35],2,mean)
mean_c_switch <- apply(moments_switch$csim[,1:35],2,mean)

mean_switch <- data.table(age = 25:59, y_switch=mean_y_switch, c_switch=mean_c_switch)

#test the real difference between rw and nl




#mean and variance of income
load('~/git/abbg/R3/rw_nbl_yparam.dat')
moments_rw <- moments

mean_y_rw <- apply(moments_rw$ysim[,1:18],2,mean)
mean_c_rw <- apply(moments_rw$csim[,1:18],2,mean)
mean_a_rw <- apply(moments_rw$asim[,1:18],2,mean)

var_y_rw <- apply(moments_rw$ysim[,1:18],2,var)
var_c_rw <- apply(moments_rw$csim[,1:18],2,var)
var_a_rw <- apply(moments_rw$asim[,1:18],2,var)

load('~/git/abbg/R3/nl_nbl_e50m50.dat')
moments_nl <- moments

mean_y_nl <- apply(moments_nl$ysim[,1:18],2,mean)
mean_c_nl <- apply(moments_nl$csim[,1:18],2,mean)
mean_a_nl <- apply(moments_nl$asim[,1:18],2,mean)

var_y_nl <- apply(moments_nl$ysim[,1:18],2,var)
var_c_nl <- apply(moments_nl$csim[,1:18],2,var)
var_a_nl <- apply(moments_nl$asim[,1:18],2,var)

mm <- data.table(age = seq(25,59,l=18),
  mean_y_rw = mean_y_rw, mean_y_nl = mean_y_nl,
  var_y_rw  = var_y_rw ,  var_y_nl  = var_y_nl ,
  mean_c_rw = mean_c_rw, mean_c_nl = mean_c_nl,
  var_c_rw  = var_c_rw ,  var_c_nl  = var_c_nl ,
  mean_a_rw = mean_a_rw, mean_a_nl = mean_a_nl,
  var_a_rw  = var_a_rw ,  var_a_nl  = var_a_nl  )


#mean and variance of eta + eps
load('~/git/abbg/R3/rw_nbl_yparam.dat')
moments_rw <- moments
moments_rw$ze <- moments_rw$zsim #+ moments_rw$esim
mean_ze_rw <- apply(moments_rw$ze[,1:18],2,mean)
var_ze_rw <- apply(moments_rw$ze[,1:18],2,var)

load('~/git/abbg/R3/nl_nbl_demean.dat')
moments_nl <- moments
moments_nl$ze <- moments_nl$zsim #+ moments_nl$esim
mean_ze_nl <- apply(moments_nl$ze[,1:18],2,mean)
var_ze_nl <- apply(moments_nl$ze[,1:18],2,var)

mm <- data.table(age = seq(25,59,l=18),
  mean_ze_rw = mean_ze_rw, mean_ze_nl = mean_ze_nl,
  var_ze_rw  = var_ze_rw ,  var_ze_nl  = var_ze_nl  )

mean_y <- mm[,1:3, with = FALSE]
mean_y_long <- reshape(mean_y, direction="long", varying=list(names(mean_y)[2:3]), v.names="Value",
        idvar=c("age"), timevar="model", times=c('rw','nl'))
p_y <- ggplot(mean_y_long, aes(x=age,y=Value)) +
          geom_line(aes(group =model, color=model)) +
          ggtitle('mean') +
          theme_bw()

ggsave('~/git/abbg/R3/figure/mean_ze.eps')

mean_y <- mm[,c(1,4:5), with = FALSE]
mean_y_long <- reshape(mean_y, direction="long", varying=list(names(mean_y)[2:3]), v.names="Value",
        idvar=c("age"), timevar="model", times=c('rw','nl'))
p_y <- ggplot(mean_y_long, aes(x=age,y=Value)) +
          geom_line(aes(group =model, color=model)) +
          ggtitle('var') +
          theme_bw()

ggsave('~/git/abbg/R3/figure/var_ze.eps')

#get eps variance from simulation of eps
rm(list=ls())
require(R.matlab)
temp <- readMat('~/git/abbg/R3/Mateps_true.mat')
Mateps <- temp[[1]]
save(Mateps, file='~/git/abbg/R3/Mateps_true.dat')

colMeans(Mateps)
veps_nl <- apply(Mateps,2,var)
#veps_nl <- lvar
save(veps_nl,file='veps_nl.dat')

mm <- data.table(age = seq(25,59,l=18),
  var_z_rw  = varzapprox ,  var_ze_nl  = varze_nl  )

veta_nl <- varzapprox
save(veta_nl,file='veta_nl.dat')