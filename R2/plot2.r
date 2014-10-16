setwd('figure')
require(ggplot2)   
require(data.table)

age = 25:94
attach(c(p,moments))

#conditional on permanent income decile, average asset
avsave_iniz <- data.table( pid = 1:nsim, iniz=zsim[,1], avsave=apply(asim[,1:70], 1, mean) )
setkey(avsave_iniz, pid)
avsave_iniz[, decile:=as.numeric(cut_number(iniz, n = 10))] #give a bin # for each person each age
setkey(avsave_iniz,decile)
avsave_iniz[,avsave_decile:=mean(avsave), by=decile]
wealth_decile <- subset(unique(avsave_iniz), select=c(decile,avsave_decile) )

#conditional on asset decile, average asset
avsave_ta <- data.table( pid = 1:nsim, ta=asim[,10], avsave=apply(asim[,10:35], 1, mean) )
setkey(avsave_ta, pid)
avsave_ta[, decile:=as.numeric(cut_number(ta, n = 10))] #give a bin # for each person each age
setkey(avsave_ta,decile)
avsave_ta[,avsave_decile:=mean(avsave), by=decile]
wealth_decile <- subset(unique(avsave_ta), select=c(decile,avsave_decile) )

######write to matlab
require(R.matlab)
load('nl.dat')
moments$lcsim <- log(moments$csim)
with(moments,writeMat('nl_zbl_risk.mat',zsim=zsim,esim=esim,asim=asim,csim=csim,lcsim=lcsim,ysim=ysim,ypresim=ypresim))

###for a given asset decile at age 34, consumption at age 35 for the 2 models? 
con_ta <- data.table( pid = 1:nsim, ta=asim[,10], con=csim[,11] )
setkey(con_ta, pid)
con_ta[, decile:=as.numeric(cut_number(ta, n = 10))] #give a bin # for each person each age
setkey(con_ta, decile)
con_ta[,con_decile:=mean(con), by=decile]
con_ta_unique <- subset(unique(con_ta), select=c(decile,con_decile) )

###for a given after tax income decile at age 34, consumption at age 35 for the 2 models? 
con_ta <- data.table( pid = 1:nsim, ta=ysim[,10], con=csim[,11] )
setkey(con_ta, pid)
con_ta[, decile:=as.numeric(cut_number(ta, n = 10))] #give a bin # for each person each age
setkey(con_ta, decile)
con_ta[,con_decile:=mean(con), by=decile]
con_ta_unique <- subset(unique(con_ta), select=c(decile,con_decile) )

###consumption between age 34 and 65, for all decile of age 34 asset
con_p <- data.table( csim[,10:41] )
con_p[,pid:=1:nsim]
con_p[,ta:=asim[,10]] 
setkey(con_p, pid)

con_p[, decile:=as.numeric(cut_number(ta, n = 10))] #give a bin # for each person each age
setkey(con_p, decile)
con_p[,pid:=NULL]
con_p[,ta:=NULL]

plot_wide <- con_p[, lapply(.SD,mean), by=decile]

plot_long <- reshape(plot_wide, direction="long", varying=list(names(plot_wide)[2:33]), v.names="Value", 
        idvar=c("decile"), timevar="age", times=34:65)
p_long <- ggplot(plot_long, aes(x=age,y=Value)) +
          geom_line(aes(group =decile, color=decile))

###consumption between age 35 and 65, for all decile of age 34 earnings
con_p <- data.table( csim[,11:41] )
con_p[,pid:=1:nsim]
con_p[,ta:=ysim[,10]] 
setkey(con_p, pid)

con_p[, decile:=as.numeric(cut_number(ta, n = 10))] #give a bin # for each person each age
setkey(con_p, decile)
con_p[,pid:=NULL]
con_p[,ta:=NULL]

plot_wide <- con_p[, lapply(.SD,mean), by=decile]

plot_long <- reshape(plot_wide, direction="long", varying=list(names(plot_wide)[2:32]), v.names="Value", 
        idvar=c("decile"), timevar="age", times=35:65)
p_long <- ggplot(plot_long, aes(x=age,y=Value)) +
          geom_line(aes(group =decile, color=decile)) +
          ggtitle('consumption conditional on earnings decile (random walk)') +
          theme_bw()
#print(p_50)
ggsave('rw.png',width=10.6, height=5.93)  

##################################################
#read matlab
require(R.matlab)
moments_nl <- readMat('nl.mat')
moments_rw <- readMat('rw.mat')
nsim <- 50000

##################################################
moments_nl <- moments
nl_fu <- with(moments_nl, data.table( pid = 1:nsim, eta=zsim[,10], con=csim[,10], ass= asim[,10] ))
rw_fu <- with(moments_rw, data.table( pid = 1:nsim, eta=zsim[,10], con=csim[,10], ass= asim[,10] ))
setkey(nl_fu, pid)
setkey(rw_fu, pid)

#eta and consumption contemporaneous 
nl_fu[, decile_fu:=as.numeric(cut_number(eta, n = 10))] #give a bin # for each person each age
rw_fu[, decile_fu:=as.numeric(cut_number(eta, n = 10))] #give a bin # for each person each age
setkey(nl_fu,  decile_fu)
setkey(rw_fu,  decile_fu)
nl_fu[,con_mean:=mean(con), by= decile_fu]
rw_fu[,con_mean:=mean(con), by= decile_fu]
nl_fu_uniq <- subset(unique(nl_fu), select=c(decile_fu,con_mean) )
rw_fu_uniq <- subset(unique(rw_fu), select=c(decile_fu,con_mean) )
nl_fu_uniq
rw_fu_uniq

nl_fu_uniq[,model:='nl']
rw_fu_uniq[,model:='rw']
fu_uniq <- rbind(nl_fu_uniq, rw_fu_uniq)

setwd('figure')
p_fu_uniq <- ggplot(fu_uniq, aes(x=decile_fu,y=con_mean)) +
          geom_line(aes(group =model, color=model)) +
          ylab('mean consumption')+
          xlab('eta decile')+
          ggtitle('consumption and eta decile at age 34') +
          theme_bw()
ggsave('con_eta_full.png',width=10.6, height=5.93)  


##################################################
#eta and consumption contemporaneous conditional on asset
nl_fu[, ass_lev:=as.numeric(cut_number(ass, n = 2))] #give a bin # for each person each age
rw_fu[, ass_lev:=as.numeric(cut_number(ass, n = 2))] #give a bin # for each person each age
nl_fu[,con_mean:=NULL]
nl_fu[,decile_fu:=NULL]
rw_fu[,con_mean:=NULL]
rw_fu[,decile_fu:=NULL]

nl_low <- nl_fu[ass_lev == 1]
nl_hig <- nl_fu[ass_lev == 2]
rw_low <- rw_fu[ass_lev == 1]
rw_hig <- rw_fu[ass_lev == 2]

nl_low[, decile:=as.numeric(cut_number(eta, n = 10))] #give a bin # for each person each age
nl_hig[, decile:=as.numeric(cut_number(eta, n = 10))] #give a bin # for each person each age
rw_low[, decile:=as.numeric(cut_number(eta, n = 10))] #give a bin # for each person each age
rw_hig[, decile:=as.numeric(cut_number(eta, n = 10))] #give a bin # for each person each age

setkey(nl_low,  decile)
setkey(nl_hig,  decile)
setkey(rw_low,  decile)
setkey(rw_hig,  decile)

nl_low[,con_mean:=mean(con), by= decile]
nl_hig[,con_mean:=mean(con), by= decile]
rw_low[,con_mean:=mean(con), by= decile]
rw_hig[,con_mean:=mean(con), by= decile]

nl_low_uniq <- subset(unique(nl_low), select=c(decile,con_mean) )
nl_hig_uniq <- subset(unique(nl_hig), select=c(decile,con_mean) )
rw_low_uniq <- subset(unique(rw_low), select=c(decile,con_mean) )
rw_hig_uniq <- subset(unique(rw_hig), select=c(decile,con_mean) )

nl_low_uniq[,model:='nl']
nl_hig_uniq[,model:='nl']
rw_low_uniq[,model:='rw']
rw_hig_uniq[,model:='rw']

nl_low_uniq[,asset:='low']
nl_hig_uniq[,asset:='high']
rw_low_uniq[,asset:='low']
rw_hig_uniq[,asset:='high']

nl_low_uniq
nl_hig_uniq
rw_low_uniq
rw_hig_uniq

uniq <- rbind(nl_low_uniq, nl_hig_uniq, rw_low_uniq, rw_hig_uniq)

p_fu_uniq <- ggplot(uniq, aes(x=decile,y=con_mean)) +
          geom_line(aes(group =interaction(model,asset), color=interaction(model,asset))) +
          ylab('mean consumption')+
          xlab('eta decile')+
          ggtitle('consumption and eta decile at age 34 conditional on asset level') +
          theme_bw()
ggsave('con_eta_asset.png',width=10.6, height=5.93)  

##################################################
moments_nl <- moments
require(R.matlab)
moments_rw <- readMat('rw.mat')

var_nl <- apply(moments_nl$ysim[,1:35],2,var)
var_rw <- apply(moments_rw$ysim[,1:35],2,var)
var_y <- data.table(age = 25:59, nl=var_nl, rw=var_rw)
var_y_long <- reshape(var_y, direction="long", varying=list(names(var_y)[2:3]), v.names="Value", 
        idvar=c("age"), timevar="model", times=c('nl','rw'))
p_y <- ggplot(var_y_long, aes(x=age,y=Value)) +
          geom_line(aes(group =model, color=model)) +
          ggtitle('variance of after tax earnings') +
          theme_bw()
#p_y
ggsave('var_earnings.png',width=10.6, height=5.93)  

mean_nl <- apply(moments_nl$ysim[,1:35],2,mean)
mean_rw <- apply(moments_rw$ysim[,1:35],2,mean)
mean_y <- data.table(age = 25:59, nl=mean_nl, rw=mean_rw)
mean_y_long <- reshape(mean_y, direction="long", varying=list(names(mean_y)[2:3]), v.names="Value", 
        idvar=c("age"), timevar="model", times=c('nl','rw'))
p_y <- ggplot(mean_y_long, aes(x=age,y=Value)) +
          geom_line(aes(group =model, color=model)) +
          ggtitle('mean of after tax earnings') +
          theme_bw()
#p_y
ggsave('mean_earnings.png',width=10.6, height=5.93)  


mm <- data.frame( age=age)
#median
mm <- cbind( mm, data.frame(stMedian = apply(asim[,-71], 2, median,na.rm=T)) )
mm <- cbind( mm, data.frame(ctMedian = apply(csim, 2, median,na.rm=T)) )
mm <- cbind( mm, data.frame(mtMedian = apply(xsim, 2, median,na.rm=T)) )
mm <- cbind( mm, data.frame(ytMedian = apply(ysim, 2, median,na.rm=T)) )

#1st quantile
mm <- cbind( mm, data.frame(st1q = apply(asim[,-71], 2, quantile, 0.25,na.rm=T)) )
mm <- cbind( mm, data.frame(ct1q = apply(csim, 2, quantile, 0.25,na.rm=T)) )
mm <- cbind( mm, data.frame(mt1q = apply(xsim, 2, quantile, 0.25,na.rm=T)) )
mm <- cbind( mm, data.frame(yt1q = apply(ysim, 2, quantile, 0.25,na.rm=T)) )

#3st quantil
mm <- cbind( mm, data.frame(st3q = apply(asim[,-71], 2, quantile, 0.75,na.rm=T)) )
mm <- cbind( mm, data.frame(ct3q = apply(csim, 2, quantile, 0.75,na.rm=T)) )
mm <- cbind( mm, data.frame(mt3q = apply(xsim, 2, quantile, 0.75,na.rm=T)) )
mm <- cbind( mm, data.frame(yt3q = apply(ysim, 2, quantile, 0.75,na.rm=T)) )

detach(moments)

mmic <- data.frame( age=age, value = mm$yt1q, moment = 'income', quantile='1st' )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ytMedian, 
  moment = 'income', quantile='Median') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$yt3q, 
  moment = 'income', quantile='3rd') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ct1q, 
  moment = 'consumption', quantile='1st') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ctMedian, 
  moment = 'consumption', quantile='Median') )
mmic <- rbind(mmic, data.frame( age=age, value = mm$ct3q, 
  moment = 'consumption', quantile='3rd') )

mma <- data.frame( age=age, asset = mm$st1q, quantile='1st') 
mma <- rbind(mma, data.frame( age=age, asset = mm$stMedian, quantile='Median') )
mma <- rbind(mma, data.frame( age=age, asset = mm$st3q, quantile='3rd') )

p_50  <- ggplot(mm, aes(x=age,y=stMedian)) + 
  geom_line(aes(color='wealth')) +
  geom_line(aes(y=ctMedian, color='consumption')) +
  #geom_line(aes(y=mtMedian, color='cash on hand')) +
  geom_line(aes(y=ytMedian, color='earnings')) +
  xlab('age') +
  ylab('value')+
  labs(colour = NULL) +
  ggtitle('life cycle medians') +
  theme_bw()
#print(p_50)
ggsave('median.png',width=10.6, height=5.93)  

p_25  <- ggplot(mm, aes(x=age,y=st1q)) + 
  geom_line(aes(color='asset')) +
  geom_line(aes(y=ct1q, color='consumption')) +
  geom_line(aes(y=mt1q, color='cash-on-hand')) +
  geom_line(aes(y=yt1q, color='income')) +
  xlab('age') +
  ylab('value')+
  labs(colour = NULL) +
  ggtitle('life cycle profile (1st quantile)') +
  theme_bw()
#print(p_25)
ggsave('1q.png',width=10.6, height=5.93)  

p_75  <- ggplot(mm, aes(x=age,y=st3q)) + 
  geom_line(aes(color='asset')) +
  geom_line(aes(y=ct3q, color='consumption')) +
  geom_line(aes(y=mt3q, color='cash-on-hand')) +
  geom_line(aes(y=yt3q, color='income')) +
  xlab('age') +
  ylab('value')+
  labs(colour = NULL) +
  ggtitle('life cycle profile (3rd quantile)') +
  theme_bw()
#print(p_75)
ggsave('3q.png',width=10.6, height=5.93)  




p_ic  <- ggplot(mmic, aes(x=age,y=value)) + 
  geom_line(aes(color=quantile,linetype=moment)) +
  xlab('age') +
  ylab('value')+
#  labs(colour = NULL) +
#  ggtitle('life type profile of income and consumption') +
  theme_bw()
print(p_ic)
ggsave('inc_cons.png',width=10.6, height=12)  

p_a  <- ggplot(mma, aes(x=age,y=asset)) + 
  geom_line(aes(color=quantile, linetype=quantile)) +
  xlab('age') +
  ylab('value')+
#  labs(colour = NULL) +
#  ggtitle('life type profile of income and consumption') +
  theme_bw()
print(p_a)
ggsave('asset.png',width=10.6, height=12) 

setwd('~/git/abbg/R')

p_inc <- ggplot(incpro,aes(x=age,y=norminc))+
  geom_line()+
  xlab('age') +
  ylab('value')+
#  labs(colour = NULL) +
  ggtitle('deterministic income profile') +
  theme_bw()
p_inc
ggsave('dinc.png',width=10.6, height=5.93)   


####plot xeta
age = seq(p$age_min, p$age_re-2, 2)
plot_xeta <- data.frame( age=age, value=xeta[,25], moment = 'node50' )
plot_xeta <- rbind( plot_xeta, data.frame( age=age, value=xeta[,1], moment = 'node1') )
plot_xeta <- rbind( plot_xeta, data.frame( age=age, value=xeta[,50], moment = 'node100') )
plot_xeta <- rbind( plot_xeta, data.frame( age=age, value=rowMeans(xeta), moment = 'means') )

p_xeta  <- ggplot(plot_xeta, aes(x=age,y=value)) + 
  geom_line(aes(color=moment)) +
  xlab('age') +
  ylab('value')+
  labs(colour = NULL) +
#  ggtitle('life type profile of income and consumption') +
  theme_bw()
print(p_xeta)
ggsave('xeta.png',width=10.6, height=12) 
