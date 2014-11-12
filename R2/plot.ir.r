require(ggplot2)   
require(data.table)
age = 25:94

load('sim.ir.res.100.0.1.0.1.dat'); ll <- moments
load('sim.ir.res.100.0.1.0.5.dat'); lm <- moments
load('sim.ir.res.100.0.1.0.9.dat'); lh <- moments
load('sim.ir.res.100.0.5.0.1.dat'); ml <- moments
load('sim.ir.res.100.0.5.0.5.dat'); mm <- moments
load('sim.ir.res.100.0.5.0.9.dat'); mh <- moments
load('sim.ir.res.100.0.9.0.1.dat'); hl <- moments
load('sim.ir.res.100.0.9.0.5.dat'); hm <- moments
load('sim.ir.res.100.0.9.0.9.dat'); hh <- moments

setwd('figure')

#pretax income
ypre <- data.table(age=age,  
  ll = colMeans(ll$ypresim),
  lm = colMeans(lm$ypresim),
  lh = colMeans(lh$ypresim), 
  ml = colMeans(ml$ypresim),
  mm = colMeans(mm$ypresim),
  mh = colMeans(mh$ypresim),
  hl = colMeans(hl$ypresim),
  hm = colMeans(hm$ypresim),
  hh = colMeans(hh$ypresim) ) 

ypre <- subset(ypre,age< 60)

ypre$dev_ll_lm <- with( ypre, (ll-lm)/lm ) 
ypre$dev_lh_lm <- with( ypre, (lh-lm)/lm ) 
ypre$dev_ml_mm <- with( ypre, (ml-mm)/mm ) 
ypre$dev_mh_mm <- with( ypre, (mh-mm)/mm ) 
ypre$dev_hl_hm <- with( ypre, (hl-hm)/hm ) 
ypre$dev_hh_hm <- with( ypre, (hh-hm)/hm ) 

p_ypre <- ggplot(ypre, aes(x=age,y=dev_ll_lm)) + 
           geom_line(aes(color='ini=.1 shock=.1')) +
           geom_line(aes(y=dev_lh_lm, color='ini=.1 shock=.9')) +
           geom_line(aes(y=dev_ml_mm, color='ini=.5 shock=.1')) +
           geom_line(aes(y=dev_mh_mm, color='ini=.5 shock=.9')) +
           geom_line(aes(y=dev_hl_hm, color='ini=.9 shock=.1')) +
           geom_line(aes(y=dev_hh_hm, color='ini=.9 shock=.9')) +
           ylab('percentage change from shock=.5')+
           labs(color = NULL) +
           ggtitle('Impulse responses, earnings') +
           theme_bw()

#p_ypre
ggsave('ir_ypre.png',width=10.6, height=5.93)  

#after tax income
y <- data.table(age=age,  
  ll = colMeans(ll$ysim),
  lm = colMeans(lm$ysim),
  lh = colMeans(lh$ysim), 
  ml = colMeans(ml$ysim),
  mm = colMeans(mm$ysim),
  mh = colMeans(mh$ysim),
  hl = colMeans(hl$ysim),
  hm = colMeans(hm$ysim),
  hh = colMeans(hh$ysim) ) 

y <- subset(y,age< 60)

y$dev_ll_lm <- with( y, (ll-lm)/lm ) 
y$dev_lh_lm <- with( y, (lh-lm)/lm ) 
y$dev_ml_mm <- with( y, (ml-mm)/mm ) 
y$dev_mh_mm <- with( y, (mh-mm)/mm ) 
y$dev_hl_hm <- with( y, (hl-hm)/hm ) 
y$dev_hh_hm <- with( y, (hh-hm)/hm ) 

p_y <- ggplot(y, aes(x=age,y=dev_ll_lm)) + 
           geom_line(aes(color='ini=.1 shock=.1')) +
           geom_line(aes(y=dev_lh_lm, color='ini=.1 shock=.9')) +
           geom_line(aes(y=dev_ml_mm, color='ini=.5 shock=.1')) +
           geom_line(aes(y=dev_mh_mm, color='ini=.5 shock=.9')) +
           geom_line(aes(y=dev_hl_hm, color='ini=.9 shock=.1')) +
           geom_line(aes(y=dev_hh_hm, color='ini=.9 shock=.9')) +
           ylab('percentage change from shock=.5')+
           labs(color = NULL) +
           ggtitle('Impulse responses, earnings') +
           theme_bw()
#p_y
ggsave('ir_ynet.png',width=10.6, height=5.93)  

#consumption
c <- data.table(age=age,  
  ll = colMeans(ll$csim),
  lm = colMeans(lm$csim),
  lh = colMeans(lh$csim), 
  ml = colMeans(ml$csim),
  mm = colMeans(mm$csim),
  mh = colMeans(mh$csim),
  hl = colMeans(hl$csim),
  hm = colMeans(hm$csim),
  hh = colMeans(hh$csim) ) 

c <- subset(c,age< 60)

c$dev_ll_lm <- with( c, (ll-lm)/lm ) 
c$dev_lh_lm <- with( c, (lh-lm)/lm ) 
c$dev_ml_mm <- with( c, (ml-mm)/mm ) 
c$dev_mh_mm <- with( c, (mh-mm)/mm ) 
c$dev_hl_hm <- with( c, (hl-hm)/hm ) 
c$dev_hh_hm <- with( c, (hh-hm)/hm ) 

p_c <- ggplot(c, aes(x=age,y=dev_ll_lm)) + 
           geom_line(aes(color='ini=.1 shock=.1')) +
           geom_line(aes(y=dev_lh_lm, color='ini=.1 shock=.9')) +
           geom_line(aes(y=dev_ml_mm, color='ini=.5 shock=.1')) +
           geom_line(aes(y=dev_mh_mm, color='ini=.5 shock=.9')) +
           geom_line(aes(y=dev_hl_hm, color='ini=.9 shock=.1')) +
           geom_line(aes(y=dev_hh_hm, color='ini=.9 shock=.9')) +
           ylab('percentage change from shock=.5')+
           labs(color = NULL) +
           ggtitle('Impulse responses, consumption') +
           theme_bw()
#p_c
ggsave('ir_con.png',width=10.6, height=5.93)  

#asset
a <- data.table(age=c(age,95),  
  ll = colMeans(ll$asim),
  lm = colMeans(lm$asim),
  lh = colMeans(lh$asim), 
  ml = colMeans(ml$asim),
  mm = colMeans(mm$asim),
  mh = colMeans(mh$asim),
  hl = colMeans(hl$asim),
  hm = colMeans(hm$asim),
  hh = colMeans(hh$asim) ) 

a <- subset(a,age>25 & age< 60)

a$dev_ll_lm <- with( a, (ll-lm)/lm ) 
a$dev_lh_lm <- with( a, (lh-lm)/lm ) 
a$dev_ml_mm <- with( a, (ml-mm)/mm ) 
a$dev_mh_mm <- with( a, (mh-mm)/mm ) 
a$dev_hl_hm <- with( a, (hl-hm)/hm ) 
a$dev_hh_hm <- with( a, (hh-hm)/hm ) 

p_a <- ggplot(a, aes(x=age,y=dev_ll_lm)) + 
           geom_line(aes(color='ini=.1 shock=.1')) +
           geom_line(aes(y=dev_lh_lm, color='ini=.1 shock=.9')) +
           geom_line(aes(y=dev_ml_mm, color='ini=.5 shock=.1')) +
           geom_line(aes(y=dev_mh_mm, color='ini=.5 shock=.9')) +
           geom_line(aes(y=dev_hl_hm, color='ini=.9 shock=.1')) +
           geom_line(aes(y=dev_hh_hm, color='ini=.9 shock=.9')) +
           ylab('percentage change from shock=.5')+
           labs(color = NULL) +
           ggtitle('Impulse responses, asset') +
           theme_bw()
#p_a
ggsave('ir_ass.png',width=10.6, height=5.93)  

