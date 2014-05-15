% Parameter

Vect=Y(:,1:T-1);
Vect_dep=Y(:,2:T);

Mat1=[];
for kk1=0:K1
    
        Mat1=[Mat1 hermite(kk1,(Vect(:)-meanY)/stdY)];
    
end

ResP_data=zeros((K1+1),Ntau);

for jtau=1:Ntau
    
    tau=Vectau(jtau);
    
    ResP_data(:,jtau)=rq(Mat1,Vect_dep(:),tau);
    
end

% Covariates (to compute the derivative of the quantile function with respect to y_{t-1})

Vect=quantile(Vect(:),Vectau);

Mat3=zeros(Ntau,1);
for kk1=1:K1
    
        Mat3=[Mat3 kk1*hermite(kk1-1,(Vect(:)-meanY)/stdY)./stdY];
   
end

% Matrix of persistence

Mat3*ResP_data

% Drawing the graph 

figure
surf(Vectau,Vectau,Mat3*ResP_data)
axis([0 1 0 1 0 1.2])
xlabel('percentile \tau_{shock}')
ylabel('percentile \tau_{init}')
zlabel('persistence')
