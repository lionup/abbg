setwd('figure')
require(ggplot2)   
age = 25:(p$NumOfPeriodsToSimulate+25-1)

plot_mean <- data.frame( cbind(age=age, 
                              income=moments$IncomeMean, 
                              consumption=moments$CtMean) )

p_mean  <- ggplot(plot_mean, aes(x=age,y=income)) + 
           geom_line(aes(color='income')) +
           geom_line(aes(y=consumption, color='consumption')) +
           xlab('age') +
           ylab('dollar')+
           labs(colour = NULL) +
           #ylab('consumption') +
           ggtitle(paste('[mean]','stdPerm=',p$sigP,'stdTran=',p$sigT)) +
           theme_bw()
ggsave('cons_inc_mean.png',width=10.6, height=5.93)           

plot_sd <- data.frame( cbind(age=age, 
                              income=apply(moments$Income,1,sd), 
                              consumption=apply(moments$Ct,1,sd)) )

p_sd  <- ggplot(plot_sd, aes(x=age,y=income)) + 
           geom_line(aes(color='income')) +
           geom_line(aes(y=consumption, color='consumption')) +
           xlab('age') +
           labs(colour = NULL) +
           #ylab('consumption') +
           ggtitle(paste('[standard deviation]','stdPerm=',p$sigP,'stdTran=',p$sigT)) +
           theme_bw()
   
#ggsave('cons_inc_sd.png',width=10.6, height=5.93)  

plot_sd_norm <- data.frame( cbind(age=age, 
          income=apply(moments$Income/moments$IncomeMean,1,sd), 
          consumption=apply(moments$Ct/moments$CtMean,1,sd)) )

p_sd_norm  <- ggplot(plot_sd_norm, aes(x=age,y=income)) + 
           geom_line(aes(color='income')) +
           geom_line(aes(y=consumption, color='consumption')) +
           xlab('age') +
           ylab('')+
           labs(colour = NULL) +
           #ylab('consumption') +
           ggtitle(paste('[normalized standard deviation]','stdPerm=',p$sigP,'stdTran=',p$sigT)) +
           theme_bw()
   
ggsave('cons_inc_sd_norm.png',width=10.6, height=5.93)            
           
setwd('~/git/abbg/R')
