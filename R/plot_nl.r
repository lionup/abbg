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




