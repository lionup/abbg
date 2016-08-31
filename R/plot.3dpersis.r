##################################################
#derivative of nonlinear persisitence
rm(list = ls())
setwd('~/git/abbg/R/figure/report14')
require(data.table)
require(EQL)
require(quantreg)
require(R.matlab)

#load('~/git/abbg/R/dat/nbin100/simdata.dat')
load('~/git/abbg/R/sim_200_199.dat')
load('~/git/abbg/R/mat.dat')

attach(data)

p <- list()
p$Resqtrue = Resqfinal
p$Resqtrue_e0 = Resqfinal.e0
p$Resqtrue_eps = Resqfinal.eps
p$b1true = b1
p$bLtrue = bL
p$b1true_e0 = b1.e0
p$bLtrue_e0 = bL.e0
p$b1true_eps = b1.eps       
p$bLtrue_eps = bL.eps

p$K1         = K1
p$K2         = K2
p$K3         = K3
p$K4         = K4
p$meanAGE    = meanAGE    
p$stdAGE     = stdAGE
p$Vectau     = Vectau
p$Ntau       = Ntau
p$meanY      = meanY
p$stdY       = stdY
p$T          = T

detach(data)

sim <- data.table(simdata)

sim$t <-1:p$T
sim_y <- sim[,c("id","t","Y"),with=F]
wide_y <- reshape(sim_y, idvar='id', timevar='t', direction='wide') #each person has all age in a row
Y <-data.matrix(wide_y[,-1,with=F])

#attach(data)
Vect=Y[,1:p$T-1]
Vect_dep=Y[,2:p$T]

#hermite decomposition of last period
Mat1<- matrix( 0, nrow=length(Vect),ncol=(p$K1+1) )
for (kk1 in 0:p$K1){
	Mat1[,kk1+1]=hermite( (c(Vect)-data$meanY)/data$stdY,kk1 )
}

#quantile regression of this period on last period hermite poly
#as hermite 0 is 1, constant, so regress without constant
ResP_data=array( 0,dim=c((p$K1+1),p$Ntau) )
Mat3 = t(ResP_data)
for (jtau in 1:p$Ntau){  
	tau=p$Vectau[jtau]
	ResP_data[,jtau]=rq(c(Vect_dep)~Mat1-1,tau)$coefficients
}

# Covariates (to compute the derivative of the quantile function with respect to y_{t-1})
# get quantiles of last period
Vect=quantile(c(Vect),p$Vectau, names=F, na.rm = T )

# hermite decompostion and derivative
# H_(n+1)(x)=xH_n(x)-n*H_(n-1)(x)
# 1, x, x^2-1, x^3-3x
# 0, 1*1, 2*x, 3*(x^2-1)
for (kk1 in 1:p$K1){
	Mat3[,kk1+1] = kk1 * hermite( (c(Vect)-data$meanY)/data$stdY,kk1-1 ) / data$stdY
}

# Matrix of persistence
persisinc <- Mat3 %*% ResP_data

names <- 'cohort30'
ename <-'_parallel' 
save( persisinc, file=paste('persisinc_',names,ename,'.dat',sep='') )
writeMat(paste('persisinc_',names,ename,'.mat',sep=''), persisinc=persisinc)


##################################################
#derivative of nonlinear persisitence
rm(list = ls())
setwd('~/git/abbg/R/figure/report14')
require(data.table)
require(EQL)
require(quantreg)
require(R.matlab)

#load('~/git/abbg/R/dat/nbin100/simdata.dat')
load('~/git/abbg/R/mat.dat')

attach(data)

p <- list()
p$Resqtrue = Resqfinal
p$Resqtrue_e0 = Resqfinal.e0
p$Resqtrue_eps = Resqfinal.eps
p$b1true = b1
p$bLtrue = bL
p$b1true_e0 = b1.e0
p$bLtrue_e0 = bL.e0
p$b1true_eps = b1.eps       
p$bLtrue_eps = bL.eps

p$K1         = K1
p$K2         = K2
p$K3         = K3
p$K4         = K4
p$meanAGE    = meanAGE    
p$stdAGE     = stdAGE
p$Vectau     = Vectau
p$Ntau       = Ntau
p$meanY      = meanY
p$stdY       = stdY
p$T          = T

detach(data)
Y <- data$Y

#attach(data)
Vect=Y[,1:p$T-1]
Vect_dep=Y[,2:p$T]

#hermite decomposition of last period
Mat1<- matrix( 0, nrow=length(Vect),ncol=(p$K1+1) )
for (kk1 in 0:p$K1){
	Mat1[,kk1+1]=hermite( (c(Vect)-data$meanY)/data$stdY,kk1 )
}

#quantile regression of this period on last period hermite poly
#as hermite 0 is 1, constant, so regress without constant
ResP_data=array( 0,dim=c((p$K1+1),p$Ntau) )
Mat3 = t(ResP_data)
for (jtau in 1:p$Ntau){  
	tau=p$Vectau[jtau]
	ResP_data[,jtau]=rq(c(Vect_dep)~Mat1-1,tau)$coefficients
}

# Covariates (to compute the derivative of the quantile function with respect to y_{t-1})
# get quantiles of last period
Vect=quantile(c(Vect),p$Vectau, names=F, na.rm = T )

# hermite decompostion and derivative
# H_(n+1)(x)=xH_n(x)-n*H_(n-1)(x)
# 1, x, x^2-1, x^3-3x
# 0, 1*1, 2*x, 3*(x^2-1)
for (kk1 in 1:p$K1){
	Mat3[,kk1+1] = kk1 * hermite( (c(Vect)-data$meanY)/data$stdY,kk1-1 ) / data$stdY
}

# Matrix of persistence
persisinc <- Mat3 %*% ResP_data

names <- 'cohort30'
ename <-'_parallel' 
save( persisinc, file=paste('persisinc_data_',names,ename,'.dat',sep='') )
writeMat(paste('persisinc_data_',names,ename,'.mat',sep=''), persisinc=persisinc)
