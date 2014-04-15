setwd('figure')
require(ggplot2)   

p_50  <- ggplot(mm, aes(x=age,y=stMedian)) + 
  geom_line(aes(color='asset')) +
  geom_line(aes(y=ctMedian, color='consumption')) +
  geom_line(aes(y=mtMedian, color='cash-on-hand')) +
  geom_line(aes(y=ytMedian, color='income')) +
  xlab('age') +
  ylab('value')+
  labs(colour = NULL) +
  ggtitle('life cycle profile (median)') +
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



##data for stephane
load('even.dat')
age = seq(p$age_min, p$age_max, 2)
attach(moments)
income <- t(ytList)
consumption <- t(ctList)
saving     <- t(stList)
eta <- t(etaList)
eps <- t(epsList)
detach(moments)
require(data.table)
sim_data <- data.table(id=1:p$nsim, age=rep(age,each=p$nsim),
            income=c(income),consumption=c(consumption),
            saving=c(saving),eta=c(eta),eps=c(eps))
setkey(sim_data,id)
sub_id <- sample(1:p$nsim,500)
sim_sub <- sim_data[J(sub_id)]
rand_age <- sample(seq(30,54,2),500,replace=T)

load('odd.dat')
age = seq(p$age_min, p$age_max, 2)
attach(moments)
income <- t(ytList)
consumption <- t(ctList)
saving     <- t(stList)
eta <- t(etaList)
eps <- t(epsList)
detach(moments)
sim_data2 <- data.table(id=(p$nsim+1):(2*p$nsim), age=rep(age,each=p$nsim),
            income=c(income),consumption=c(consumption),
            saving=c(saving),eta=c(eta),eps=c(eps))

sim <- rbind(sim_data,sim_data2)
setkey(sim,id,age)
sub_id <- sample(1:(2*p$nsim),1000)
sim_sub <- sim[J(sub_id)]

#simdata <-data.matrix(sim_data)
#save(simdata,file='simdata.dat')
#writeMat('simdata.mat',simdata=simdata)