ThreewayME.f <- function(M,X,Z,W,xlab,ylab,lloc,Min,Q1,Mean,Q3,Max,level,rob,hist)  
{ 
  ## A function to generate 3-way marginal effects plots in R with stars for set confidence levels. 
  ## Written by Joshua Gubler ~ http://scholar.byu.edu/jgubler
  ## Last modified: 25 June 2014
  ## Variables must be in the following order: y = x z w (control variables here) xz xw zw xzw .  Of course, the model can include as many control variables as you need; the following code specifically uses just the variables we want for the interaction.  As long as you get the first three variables in correct order, R will correctly order the interaction terms.
  ## M = an object of type "lm," "glm," or other estimation -- i.e. the object that contains the regression estimation you seek to plot
  ## X = the variable whose effect on Y you seek to plot
  ## Z = the first moderating variable (will be positioned on the X-axis of the plot)
  ## W = the second moderating variable (the lines on the plot)
  ## xlab = Label for x-axis (in quotes)
  ## ylab = label for y-axis (in quotes)
  ## lloc = location of the legend for the plot, use values like "bottomleft"
  ## Min, Q1, Mean, Q3, Max = titles for each of these quartiles to be put in the legend -- i.e. "Min Q88" (titles must be in quotes)
  ## level = to set the confidence level.  Two options (don't put these in quotes): 95, 90.  Stars will show on lines that are significant at the level you set.  If you do not put in either option, stars will show on all lines.
  ## rob = if you desire Huber-White robust standard errors
  ## hist = if hist, then you will get a density histogram of your variable on the X-axis.  If not, then you will get a rug plot of your variable on the X-axis
  
  ## Example: ThreewayME.f(estimation.lm,ses,edu,pop,"Education levels","Effect of SES on Civil War","bottomleft","Min Pop","1Q Pop","Mean Pop","3Q Pop","Max Pop",90,rob,hist)

  # n <- 25
  # x <- rexp(n)
  # z <- rexp(n)
  # w <- rexp(n)                                                                  
  # y <- x + z + x*z + x*w + rnorm(n)
  # m <- lm(y ~ x + z + x*z + x*w)
  
  S <- summary(M)
  N <- c(1:20)
  
  #     ################################################################  #  
  #       Create 20 equally spaced values in a vector between min         #
  #       and max on the Z variable                                       # 
  #     ################################################################  #
  
  zmin <- rep(min(Z,na.rm=TRUE), 20) 
  zmax <- rep(max(Z, na.rm=TRUE), 20) 
  Znew <- (((N-1)/(20-1))*(zmax-zmin))+zmin
  
  #     ################################################################  # 
  #       Generate the values of W for which you want to calculate the    #
  #       marginal effect (and standard errors) of X on Y.                #
  #     ################################################################  #
  
  W0 <- quantile(W,   0, na.rm=TRUE) 
  W1 <- quantile(W, .25, na.rm=TRUE)
  W2 <- quantile(W, .50, na.rm=TRUE)
  W3 <- quantile(W, .75, na.rm=TRUE)
  W4 <- quantile(W,   1, na.rm=TRUE)
  
  #     ################################################################  #
  #       Grab elements of the coefficient and variance-covariance matrix #
  #       that are required to calculate the marginal effect and standard #
  #       errors.                                                         #
  #     ################################################################  #
  
  require(sandwich)
  H <- head(S$coefficients,4)
  T <- tail(S$coefficients,4)
  b <- rbind(H,T)
  if(missing(rob)){Vcov <- vcov(M)}
  else{
    X <- model.matrix(M)
    u2 <- residuals(M)^2
    XDX <- 0
    for(i in 1:nrow(X)) {
      XDX <- XDX + u2[i]*X[i,]%*%t(X[i,])
    }
    XX1 <- solve(t(X)%*%X)
    Vcov <- XX1 %*% XDX %*% XX1
    dfc <- sqrt(nrow(X))/sqrt(nrow(X)-ncol(X))
  }
  Vcov <- as.data.frame(Vcov)
  Vcov1 <- Vcov[,c(1:4)]
  Vcov2 <- Vcov[,-c(3:0-length(Vcov))]
  Vcov <- cbind(Vcov1,Vcov2)
  Vh <- head(Vcov,4)
  Vt <- tail(Vcov,4)
  V <- rbind(Vh,Vt)
  
  b1 <- b[2,1]
  b2 <- b[3,1]
  b3 <- b[4,1]
  b4 <- b[5,1]
  b5 <- b[6,1]
  b6 <- b[7,1]
  b7 <- b[8,1]
  
  varb1 <- V[2,2]
  varb2 <- V[3,3]
  varb3 <- V[4,4]
  varb4 <- V[5,5]
  varb5 <- V[6,6]
  varb6 <- V[7,7]
  varb7 <- V[8,8]
  
  covb1b4 <- V[5,2]
  covb1b5 <- V[6,2]
  covb1b7 <- V[8,2]
  covb4b5 <- V[6,5]
  covb4b7 <- V[8,5]
  covb5b7 <- V[8,6]
  
  #     ################################################################  #
  #       We want to calculate the marginal effect of X on Y for all      #
  #       Z values of the modifying variable Z. We also want to           #
  #       calculate this marginal effect as Z changes for specific values #
  #       of the second modifying variable W.  In the code below, we      #
  #       calculate the marginal effect of X on Y for all values of Z     #
  #       when W=0, when W=1, when W=2, when W=3, and when W=4.           #
  #     ################################################################  #
  
  conb0 <- b1+b4*Znew+b5*W0+b7*(Znew*W0)
  conb1 <- b1+b4*Znew+b5*W1+b7*(Znew*W1)
  conb2 <- b1+b4*Znew+b5*W2+b7*(Znew*W2)
  conb3 <- b1+b4*Znew+b5*W3+b7*(Znew*W3)
  conb4 <- b1+b4*Znew+b5*W4+b7*(Znew*W4)
  
  #     ################################################################  #
  #       Calculate the standard errors for the marginal effect of X on Y #
  #       for all Z values of the modifying variable Z. Do this for the   #
  #       case when W=0, when W=1, when W=2, when W=3, and when W=4.      #
  #     ################################################################  #
  
  if(missing(rob)){
  conse0 <- sqrt(varb1         
                 + varb4*(Znew^2) + varb5*(W0^2) + varb7*(Znew^2)*(W0^2)
                 + 2*Znew*covb1b4 + 2*W0*covb1b5 + 2*Znew*W0*covb1b7 + 2*Znew*W0*covb4b5
                 + 2*W0*(Znew^2)*covb4b7 + 2*(W0^2)*Znew*covb5b7)
  
  conse1 <- sqrt(varb1         
                 + varb4*(Znew^2) + varb5*(W1^2) + varb7*(Znew^2)*(W1^2)      
                 + 2*Znew*covb1b4 + 2*W1*covb1b5 + 2*Znew*W1*covb1b7 + 2*Znew*W1*covb4b5
                 + 2*W1*(Znew^2)*covb4b7 + 2*(W1^2)*Znew*covb5b7)
  
  conse2 <- sqrt(varb1       
                 + varb4*(Znew^2) + varb5*(W2^2) + varb7*(Znew^2)*(W2^2) 
                 + 2*Znew*covb1b4 + 2*W2*covb1b5 + 2*Znew*W2*covb1b7 + 2*Znew*W2*covb4b5
                 + 2*W2*(Znew^2)*covb4b7 + 2*(W2^2)*Znew*covb5b7)
  
  conse3 <- sqrt(varb1
                 + varb4*(Znew^2) + varb5*(W3^2) + varb7*(Znew^2)*(W3^2)
                 + 2*Znew*covb1b4 + 2*W3*covb1b5 + 2*Znew*W3*covb1b7 + 2*Znew*W3*covb4b5
                 + 2*W3*(Znew^2)*covb4b7 + 2*(W3^2)*Znew*covb5b7)      
  
  conse4 <- sqrt(varb1
                 + varb4*(Znew^2) + varb5*(W4^2) + varb7*(Znew^2)*(W4^2)
                 + 2*Znew*covb1b4 + 2*W4*covb1b5 + 2*Znew*W4*covb1b7 + 2*Znew*W4*covb4b5
                 + 2*W4*(Znew^2)*covb4b7 + 2*(W4^2)*Znew*covb5b7)    
  }else{
    conse0 <- dfc*sqrt(varb1         
                   + varb4*(Znew^2) + varb5*(W0^2) + varb7*(Znew^2)*(W0^2)
                   + 2*Znew*covb1b4 + 2*W0*covb1b5 + 2*Znew*W0*covb1b7 + 2*Znew*W0*covb4b5
                   + 2*W0*(Znew^2)*covb4b7 + 2*(W0^2)*Znew*covb5b7)
    
    conse1 <- dfc*sqrt(varb1         
                   + varb4*(Znew^2) + varb5*(W1^2) + varb7*(Znew^2)*(W1^2)      
                   + 2*Znew*covb1b4 + 2*W1*covb1b5 + 2*Znew*W1*covb1b7 + 2*Znew*W1*covb4b5
                   + 2*W1*(Znew^2)*covb4b7 + 2*(W1^2)*Znew*covb5b7)
    
    conse2 <- dfc*sqrt(varb1       
                   + varb4*(Znew^2) + varb5*(W2^2) + varb7*(Znew^2)*(W2^2) 
                   + 2*Znew*covb1b4 + 2*W2*covb1b5 + 2*Znew*W2*covb1b7 + 2*Znew*W2*covb4b5
                   + 2*W2*(Znew^2)*covb4b7 + 2*(W2^2)*Znew*covb5b7)
    
    conse3 <- dfc*sqrt(varb1
                   + varb4*(Znew^2) + varb5*(W3^2) + varb7*(Znew^2)*(W3^2)
                   + 2*Znew*covb1b4 + 2*W3*covb1b5 + 2*Znew*W3*covb1b7 + 2*Znew*W3*covb4b5
                   + 2*W3*(Znew^2)*covb4b7 + 2*(W3^2)*Znew*covb5b7)      
    
    conse4 <- dfc*sqrt(varb1
                   + varb4*(Znew^2) + varb5*(W4^2) + varb7*(Znew^2)*(W4^2)
                   + 2*Znew*covb1b4 + 2*W4*covb1b5 + 2*Znew*W4*covb1b7 + 2*Znew*W4*covb4b5
                   + 2*W4*(Znew^2)*covb4b7 + 2*(W4^2)*Znew*covb5b7)    
  }
   
  #     ################################################################  #
  #                           Create t statistics                         #
  #     ################################################################  #
  
  t0 <- conb0/conse0
  t1 <- conb1/conse1
  t2 <- conb2/conse2
  t3 <- conb3/conse3
  t4 <- conb4/conse4
  
  #     ################################################################  #
  #       Make a `shadow' variable that is missing if the t score is not  #
  #       larger than the critical level of significance you want.        #
  #     ################################################################  #
  
  ci <- NA
  ci[level==95] <- 1.96
  ci[level==90] <- 1.645
  
  stars.df <- data.frame(consb0=conb0,consb1=conb1,consb2=conb2,consb3=conb3,consb4=conb4,t0=t0,t1=t1,t2=t2,t3=t3,t4=t4)
  
  stars.df$consb0[abs(stars.df$t0)<ci] <- NA
  stars.df$consb1[abs(stars.df$t1)<ci] <- NA
  stars.df$consb2[abs(stars.df$t2)<ci] <- NA
  stars.df$consb3[abs(stars.df$t3)<ci] <- NA
  stars.df$consb4[abs(stars.df$t4)<ci] <- NA
  
  #     ################################################################  #
  #       Generate a string variable called txt that is designated with  #
  #       a star.                                                         #
  #     ################################################################  #
  
  txt <- c("*")
  
  #     ################################################################  #
  #       Graph the marginal effect of X on Y across the desired range of #
  #       the modifying variable Z. Do this for when W=0, when W=1, when  #
  #       W=2, when W=3, and when W=4.                                    #
  #     ################################################################  #
  
  if(missing(hist)){
  plot(c(Znew,Znew,Znew,Znew,Znew,Znew,Znew,Znew,Znew,Znew), c(conb0,stars.df$consb0,conb1,stars.df$consb1,conb2,stars.df$consb2,conb3,stars.df$consb3,conb4,stars.df$consb4), type="n",xlab=xlab,ylab=ylab)+rug(jitter(Z))
  
  lines(Znew,conb0,col="blue")
  text(x=Znew,y=stars.df$consb0,labels=txt)
  lines(Znew,conb1,col="red")
  text(x=Znew,y=stars.df$consb1,labels=txt)
  lines(Znew,conb2,col="forest green")
  text(x=Znew,y=stars.df$consb2,labels=txt)
  lines(Znew,conb3,col="yellow")
  text(x=Znew,y=stars.df$consb3,labels=txt)
  lines(Znew,conb4,col="brown")
  text(x=Znew,y=stars.df$consb4,labels=txt)
  legend(lloc,legend=c(Min,Q1,Mean,Q3,Max),col=c("blue","red","forest green","yellow","brown"),lty = c("solid"))
  abline(h=0)
  }else{
  hist(Z, axes=F, xlab="", ylab="",main="", col="light gray")
  par(new=TRUE)
  plot(c(Znew,Znew,Znew,Znew,Znew,Znew,Znew,Znew,Znew,Znew), c(conb0,stars.df$consb0,conb1,stars.df$consb1,conb2,stars.df$consb2,conb3,stars.df$consb3,conb4,stars.df$consb4), type="n",xlab=xlab,ylab=ylab)
  
  lines(Znew,conb0,col="blue")
  text(x=Znew,y=stars.df$consb0,labels=txt)
  lines(Znew,conb1,col="red")
  text(x=Znew,y=stars.df$consb1,labels=txt)
  lines(Znew,conb2,col="forest green")
  text(x=Znew,y=stars.df$consb2,labels=txt)
  lines(Znew,conb3,col="yellow")
  text(x=Znew,y=stars.df$consb3,labels=txt)
  lines(Znew,conb4,col="brown")
  text(x=Znew,y=stars.df$consb4,labels=txt)
  legend(lloc,legend=c(Min,Q1,Mean,Q3,Max),col=c("blue","red","forest green","yellow","brown"),lty = c("solid"))
  abline(h=0)
  }
}
