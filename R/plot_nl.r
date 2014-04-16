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

p_inc <- ggplot(incpro,aes(x=age,y=norminc))+
  geom_line()+
  xlab('age') +
  ylab('value')+
#  labs(colour = NULL) +
  ggtitle('deterministic income profile') +
  theme_bw()
p_inc
ggsave('dinc.png',width=10.6, height=5.93)   

##data for stephane
load('even.dat')
age_full = seq(p$age_min, p$age_max, 2)
last_age = p$age_re - 2*6
age_nore = seq(p$age_min, last_age, 2)
attach(moments)
income <- t(ytList)
consumption <- t(ctList)
saving     <- t(stList)
eta <- t(etaList)
eps <- t(epsList)
detach(moments)
require(data.table)
sim_data <- data.table(id=1:p$nsim, age=rep(age_full,each=p$nsim),
            income=c(income),consumption=c(consumption),
            saving=c(saving),eta=c(eta),eps=c(eps))
setkey(sim_data,id)

#for each age, randomly pick 40 people
sub_id <- array( sample(1:p$nsim, length(age_nore)*40), 
  dim=c(length(age_nore),40) )

sim <- data.table()
for( t in 1:length(age_nore) ){
  sim_sub <- sim_data[J(sub_id[t,])]
  sim_sub <- subset(sim_sub, age>=age_nore[t] & age<=age_nore[t]+10)
  sim <- rbind(sim,sim_sub)
}
#######################################
load('odd.dat')
age_full = seq(p$age_min, p$age_max, 2)
last_age = p$age_re - 2*6
age_nore = seq(p$age_min, last_age, 2)
attach(moments)
income <- t(ytList)
consumption <- t(ctList)
saving     <- t(stList)
eta <- t(etaList)
eps <- t(epsList)
detach(moments)
sim_data2 <- data.table(id=(p$nsim+1):(2*p$nsim), age=rep(age_full,each=p$nsim),
            income=c(income),consumption=c(consumption),
            saving=c(saving),eta=c(eta),eps=c(eps))
setkey(sim_data2,id)

#for each age, randomly pick 40 people
sub_id <- array( sample((p$nsim+1):(2*p$nsim), length(age_nore)*40), 
  dim=c(length(age_nore),40) )

for( t in 1:length(age_nore) ){
  sim_sub <- sim_data2[J(sub_id[t,])]
  sim_sub <- subset(sim_sub, age>=age_nore[t] & age<=age_nore[t]+10)
  sim <- rbind(sim,sim_sub)
}
setkey(sim,id,age)

simdata <-data.matrix(sim)
save(simdata,file='simdata.dat')
require(R.matlab)
writeMat('simdata.mat',simdata=simdata)