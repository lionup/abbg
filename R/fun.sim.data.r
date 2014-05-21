##data for stephane
sim.origin.sample <- function(){
  load('even.dat')
  age_full = seq(p$age_min, p$age_max, 2)
  last_age = p$age_re - 2*p$T  #oldest cohort ini age
  age_nore = seq(p$age_min, last_age, 2)

  require(data.table)
  sim_data <- with( moments, data.table(id=rep(1:p$nsim,each=p$nage+p$ntr), 
    age=age_full, income=c(ytList),consumption=c(ctList),
    saving=c(stList),eta=c(etaList),eps=c(epsList)) )
  sim_data[,Y:=eta+eps]
  setkey(sim_data,id)

  #for each age, randomly pick people
  nsim_age <- trunc( p$nsim/length(age_nore) )
  sub_id <- array( sample(1:p$nsim, nsim_age*length(age_nore)), 
    dim=c(length(age_nore),nsim_age) )

  sim <- data.table()
  for( t in 1:length(age_nore) ){
    sim_sub <- sim_data[J(sub_id[t,])]
    sim_sub <- subset( sim_sub, age>=age_nore[t] & age<=age_nore[t]+2*c(p$T)-2 )
    sim <- rbind(sim,sim_sub)
  }
  #######################################
  load('odd.dat')
  age_full = seq(p$age_min, p$age_max, 2)
  last_age = p$age_re - 2*p$T  #oldest cohort ini age
  age_nore = seq(p$age_min, last_age, 2)

  sim_data2 <- with( moments, data.table(id=rep((p$nsim+1):(2*p$nsim),each=p$nage+p$ntr), 
    age=age_full, income=c(ytList),consumption=c(ctList),
    saving=c(stList),eta=c(etaList),eps=c(epsList)) )
  sim_data2[,Y:=eta+eps]
  setkey(sim_data2,id)

  #for each age, randomly pick 40 people
  nsim_age <- trunc( p$nsim/length(age_nore) )
  sub_id <- array( sample((p$nsim+1):(2*p$nsim), nsim_age*length(age_nore)), 
    dim=c(length(age_nore),nsim_age) )

  for( t in 1:length(age_nore) ){
    sim_sub <- sim_data2[J(sub_id[t,])]
    sim_sub <- subset(sim_sub, age>=age_nore[t] & age<=age_nore[t]+2*c(p$T)-2 )
    sim <- rbind(sim,sim_sub)
  }
  setkey(sim,id,age) 
  return(sim)
}


sim.small.sample <- function(iniage){
  savename <- paste('cohort',iniage,'.dat',sep='')
  load(savename)
  age_full = seq(p$age_min, p$age_max, 2)

  require(data.table)
  sim_data <- with( etaeps, data.table(id=rep(1:p$nsim,each=p$nage), 
    age=age_full, eta=c(etaList),eps=c(epsList)) )
  sim_data[,Y:=eta+eps]
  setkey(sim_data,id)
  return(sim_data)
}



#################################################
#plot persistent 
###################################################
sim.persis <- function(p, sim){

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
  require(quantreg)
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
  persis <- Mat3 %*% ResP_data

  # Drawing the graph 
  require(plot3D)
  #png('figure/persis_y.png',width=10.6, height=5.93, units='in', res=300)
  persp3D(x=c(p$Vectau),y=c(p$Vectau),z=persis, 
    xlab='percentile initial', ylab='percentile shock', zlab='persistence', 
    ticktype = "detailed")
  #dev.off() 
}