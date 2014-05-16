sim.qr <- function(p, sim){

sim$t <-1:6 
sim_y <- sim[,c("id","t","Y"),with=F]
wide_y <- reshape(sim_y, idvar='id', timevar='t', direction='wide') #each person has all age in a row

Vect=Y[,1:T-1]
Vect_dep=Y[,2:T]

Mat1<- matrix( 0, nrow=length(Vect),ncol=(K1+1) )
for (kk1 in 0:K1){
  Mat1[,kk1+1]=hermite( (c(Vect)-mean(Y))/sd(Y),kk1 )
}

ResP_data=array( 0,dim=c((K1+1),Ntau) )

require(quantreg)
for (jtau in 1:Ntau){
    
    tau=Vectau[jtau]
    
    ResP_data[,jtau]=rq(c(Vect_dep)~Mat1-1,tau)$coefficients
    
}

# Covariates (to compute the derivative of the quantile function with respect to y_{t-1})

Vect=quantile(c(Vect),Vectau, names=F, na.rm = T )

Mat3=matrix( 1, nrow=Ntau,ncol=K1+1 )
for (kk1 in 1:K1){
  Mat3[,kk1+1] = kk1 * hermite( (c(Vect)-mean(Y))/sd(Y),kk1-1 ) / sd(Y)
}

# Matrix of persistence

persis <- Mat3 %*% ResP_data

# Drawing the graph 

require(plot3D)
persp3D(Vectau,Vectau,persis)
axis([0 1 0 1 0 1.2])
xlabel('percentile \tau_{shock}')
ylabel('percentile \tau_{init}')
zlabel('persistence')

}