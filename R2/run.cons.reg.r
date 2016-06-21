##################################################
#derivative of consumption on earnings, by asset and by age
rm(list = ls())
setwd('~/git/abbg/R2/figure/report14')
require(ggplot2)   
require(data.table)
require(EQL)

names <- 'nl_nbl'
#ename <- '_etaeps'
ename <- ''
load( paste('~/git/abbg/R2/',names,'.dat',sep='') )
moments$lcsim <- log(moments$csim)
#moments$asim[moments$asim<1e-12] <- 1e-12
#moments$lasim <- log(moments$asim)
moments$lysim <- log(moments$ysim)

#Q1: is working age correct? 
attach(c(p,moments))
age = 25:(25+Twork-1)
nobs = nsim*Twork
nl_fu <- with(moments, data.table( pid = 1:nsim,  age=rep(age,each=nsim), 
	eta=c(zsim), eps=c(esim), inc=c(lysim[1:nobs]), con=c(lcsim[1:nobs]), 
	ass=c(asim[1:nobs]) ))
nl_fu <- nl_fu[ass>0] 
nl_fu[,ass:=log(ass)]
setkey(nl_fu, pid, age)

#take out full age dummies for consumption and asset
nl_fu[ ,cres:=lm(con~factor(age))$residuals ]
nl_fu[ ,ares:=lm(ass~factor(age))$residuals ]
nl_fu[ ,yres:=lm(inc~factor(age))$residuals ]
#nl_fu[,yres:=eta+eps]


#cut data into asset and age deciles
#nl_fu[, agedec:=as.numeric(cut_number(age, n = 10))] #give a bin # for each person each age
#nl_fu[, assdec:=as.numeric(cut_number(ares, n = 10))] #give a bin # for each person each age

K1 <- 2 #second order hermite for income
K2 <- 1 #second order hermite for age
K3 <- 2 #second order hermite for ass
ntau <- 11  #number of tau used as interpolation node
VecTau <- (1:ntau)/(1+ntau) #get the quantile of interpolation node

#standarlize consumption, income, age, and 
#nl_fu[,cstd:=(cres-mean(cres))/sd(cres)]
nl_fu[,astd:=(ares-mean(ares))/sd(ares)]
nl_fu[,ystd:=(yres-mean(yres))/sd(yres)]
nl_fu[,agestd:=(age-mean(age))/sd(age)]

#evaluate c=g(y,a,age)
Mat = array( 0, dim=c(nrow(nl_fu),(K1+1)*(K2+1)*(K3+1)) )

for (kk1 in 0:K1) {
  for (kk2 in 0:K2) {    
    for (kk3 in 0:K3) {
      Mat[,1+kk3+(K3+1)*(kk2+(K2+1)*kk1)] = hermite(nl_fu$ystd, kk1) * 
      	hermite(nl_fu$agestd, kk2) * hermite(nl_fu$astd, kk3)
    }
  }  
}
save( Mat, file=paste('Mat_',names,ename,'.dat',sep='') )
#load( paste('Mat_',names,ename,'.dat',sep='') )

#require(parallel)
#decomHerm <- function(i, K1, K2, K3, nl_fu){
#	for (kk1 in 0:K1) {
#	  for (kk2 in 0:K2) {    
#	    for (kk3 in 0:K3) {
#	      Mat[i,1+kk3+(K3+1)*(kk2+(K2+1)*kk1)] = hermite(nl_fu$ystd[i], kk1) * 
#	      	hermite(nl_fu$agestd[i], kk2) * hermite(nl_fu$astd[i], kk3)
#	    }
#	  }  
#	}
#}
#Mat <- mclapply(1:nrow(nl_fu), decomHerm, K1, K2, K3, nl_fu)

#regress cons on inc ass age 
#Q2: do I have to do fixed effects? 
ResP <- lm(cres~Mat-1,data=nl_fu)$coefficients
#ResP <- lm(cstd~Mat-1,data=nl_fu)$coefficients

# Covariates (to compute the derivative of the quantile function with respect to y_{t-1})
# get quantiles of last period
Vage = quantile(c(nl_fu$agestd),VecTau, names=F, na.rm = T ) 
Vass = quantile(c(nl_fu$astd),VecTau, names=F, na.rm = T )
meaninc = mean(nl_fu$ystd) #this equal to 0
sdinc = sd(nl_fu$yres)
Vgrid = data.matrix( expand.grid(z1=Vage, z2=Vass) )

# hermite decompostion and derivative
# H_(n+1)(x)=xH_n(x)-n*H_(n-1)(x)
# 1, x, x^2-1, x^3-3x
# 0, 1*1, 2*x, 3*(x^2-1)
Mat3 = array( 0, dim=c(ntau^2, (K1+1)*(K2+1)*(K3+1)) )
for (kk1 in 1:K1) {
  for (kk2 in 0:K2) {    
    for (kk3 in 0:K3) {
		  Mat3[,1+kk3+(K3+1)*(kk2+(K2+1)*kk1)] = kk1 * hermite( meaninc,kk1-1 ) / sdinc* 
		    hermite(Vgrid[,1], kk2) * hermite(Vgrid[,2], kk3)
    }
  }  
}

# Matrix of persistence
persis2 <- Mat3 %*% ResP
persis <- array(persis, dim=c(ntau, ntau))

require(plot3D)
png(paste('dc_',names,ename,'.png',sep=''),width=10.6, height=5.93, units='in', res=300)

persp3D(x=VecTau,y=VecTau,z=persis, 
  xlab='age percentile', ylab='asset percentile', zlab='consumption response', zlim=c(0,0.8),
  ticktype = "detailed",theta=-60, phi=25)

dev.off() 


setkey(nl_fu,assdec, agedec)
nl_fu[,dcdi:=lm(cres~inc)$coefficients[2],by=list(assdec, agedec)]

#keep unique number
wealth_decile <- subset(unique(nl_fu), select=c(assdec,agedec,dcdi) )

