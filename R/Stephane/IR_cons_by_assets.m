
clear all
clc;

close all


load('data_hermite_cons.mat');

Resqtrue=Resqfinal;

Resqtrue_e0=Resqfinal_e0;

Resqtrue_eps=Resqfinal_eps;

b1true=b1;
bLtrue=bL;
b1true_e0=b1_e0;
bLtrue_e0=bL_e0;
b1true_eps=b1_eps;
bLtrue_eps=bL_eps;

Restrue=Resfinal;

Restrue_a1=Resfinal_a1;

b1true_a=b1_a;
bLtrue_a=bL_a;


% seeds
rng('shuffle')

% Number of individuals
N=100000;

% age
% aa_ref=47;
% 
% nage=9;

%aa_ref=35;

% nage=15;

aa_ref=53;

nage=6;



% Keep only tau1-percentile of initial eta

tau1=.10;

% Proposal, eta_1
V_draw=unifrnd(0,1,N,1);

MatAGE1=[];
for kk3=0:K3
    MatAGE1=[MatAGE1 hermite(kk3,(aa_ref-meanAGE)/stdAGE)];
end

Mateta_true=zeros(N,nage);
Mateps_true=zeros(N,nage);
Y=zeros(N,nage);
C=zeros(N,nage);
A=zeros(N,nage);



%First quantile
Mateta_true(:,1)=(MatAGE1*Resqtrue_e0(:,1)).*(V_draw<=Vectau(1));
for jtau=2:Ntau
    Mateta_true(:,1)=Mateta_true(:,1)+((MatAGE1*(Resqtrue_e0(:,jtau)-Resqtrue_e0(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
        (V_draw-Vectau(jtau-1))+MatAGE1*Resqtrue_e0(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
end
%Last quantile.
Mateta_true(:,1)=Mateta_true(:,1)+(MatAGE1*Resqtrue_e0(:,Ntau)).*(V_draw>Vectau(Ntau));

% Laplace tails
Mateta_true(:,1)=Mateta_true(:,1)+((1/(b1true_e0)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
    -(1/bLtrue_e0*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));

% Initial eta

eta_1=quantile(Mateta_true(:,1),tau1);

Mateta_true(:,1)=eta_1;

% Initial assets


Mateta1=[];
for mm5=0:M5
    for mm6=0:M6
        Mateta1=[Mateta1 hermite(mm5,(Mateta_true(:,1)-meanY)/stdY).*hermite(mm6,(aa_ref-meanAGE)/stdAGE)];
    end
end




% Percentile initial assets
tau3=.10;
V_draw=tau3;

% Proposal, a1


%First quantile
A(:,1)=(Mateta1*Restrue_a1(:,1)).*(V_draw<=Vectau(1));
for jtau=2:Ntau
    A(:,1)=A(:,1)+((Mateta1*(Restrue_a1(:,jtau)-Restrue_a1(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
        (V_draw-Vectau(jtau-1))+Mateta1*Restrue_a1(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
end
%Last quantile.
A(:,1)=A(:,1)+(Mateta1*Restrue_a1(:,Ntau)).*(V_draw>Vectau(Ntau));

% Laplace tails
A(:,1)=A(:,1)+((1/(b1true_a)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
    -(1/bLtrue_a*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));




%%%% FIRST SIMULATION

% percentile of initial shock
tau0=.50;

for jj=1:nage
    
    aa=aa_ref+(jj-1)*2;
    
    MatAGE1=[];
    for kk3=0:K3
        MatAGE1=[MatAGE1 hermite(kk3,(aa-meanAGE)/stdAGE)];
    end
    
    % epsilons
    
    MatAGE_t=[];
    for kk4=0:K4
        MatAGE_t=[MatAGE_t hermite(kk4,(aa-meanAGE)/stdAGE)];
    end
    
    
    
    % Proposal, eta_0
    V_draw=unifrnd(0,1,N,1);
    
    %First quantile
    Mateps_true(:,jj)=(MatAGE_t*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
    for jtau=2:Ntau
        Mateps_true(:,jj)=Mateps_true(:,jj)+((MatAGE_t*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+MatAGE_t*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    end
    %Last quantile.
    Mateps_true(:,jj)=Mateps_true(:,jj)+(MatAGE_t*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
    
    % Laplace tails
    Mateps_true(:,jj)=Mateps_true(:,jj)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
    
    
    % Earnings
    
    Y(:,jj)=Mateta_true(:,jj)+Mateps_true(:,jj);
    
    
    
    
    % Assets and consumption
    
    
    XX=[];
    for mm1=0:M1
        for mm2=0:M2
            for mm3=0:M3
                for mm4=0:M4
                    XX=[XX hermite(mm1,(A(:,jj)-meanA)/stdA).*hermite(mm2,(Mateta_true(:,jj)-meanY)/stdY).*hermite(mm3,(Y(:,jj)-Mateta_true(:,jj)-meanY)/stdY).*hermite(mm4,(aa-meanAGE)/stdAGE)];
                end
            end
        end
    end
    
    C(:,jj)=XX*Restrue(1:(M1+1)*(M2+1)*(M3+1)*(M4+1),1)+sqrt(Restrue((M1+1)*(M2+1)*(M3+1)*(M4+1)+1))*randn(N,1);
    
    
    
    if jj<=nage-1
        
        
        %A(:,jj+1)=[ones(N,1) A(:,jj) Y(:,jj) C(:,jj)]*Res_assets+sqrt(sig_assets)*randn(N,1);
        A(:,jj+1)=log(max((1+.03)*exp(A(:,jj))+exp(Y(:,jj))-exp(C(:,jj)),ones(N,1)));
        
        Mat=zeros(N,(K1+1)*(K2+1));
        for kk1=0:K1
            for kk2=0:K2
                Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true(:,jj)-meanY)/stdY).*hermite(kk2,((aa+2)-meanAGE)/stdAGE);
            end
        end
        
        % age 37 receives shock tau0
        if jj==1
            V_draw=tau0*ones(N,1);
        else
            V_draw=unifrnd(0,1,N,1);
        end
        
        %First quantile
        
        Mateta_true(:,jj+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1));
        for jtau=2:Ntau
            Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+...
                ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
                (Vectau(jtau)-Vectau(jtau-1)).*...
                (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
                (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
        end
        %Last quantile.
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+(Mat*Resqtrue(:,Ntau)).*...
            (V_draw>Vectau(Ntau));
        
        % Laplace tails
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
            -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
        
    end
end

Y1=Y;
C1=C;
Mateta_true1=Mateta_true;
A1=A;

%%%% SECOND SIMULATION

% percentile of initial shock
tau0=.10;

for jj=1:nage
    
    aa=aa_ref+(jj-1)*2;
    
   
    
    MatAGE1=[];
    for kk3=0:K3
        MatAGE1=[MatAGE1 hermite(kk3,(aa-meanAGE)/stdAGE)];
    end
    
    % epsilons
    
    MatAGE_t=[];
    for kk4=0:K4
        MatAGE_t=[MatAGE_t hermite(kk4,(aa-meanAGE)/stdAGE)];
    end
    
    
    
    % Proposal, eta_0
    V_draw=unifrnd(0,1,N,1);
    
    %First quantile
    Mateps_true(:,jj)=(MatAGE_t*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
    for jtau=2:Ntau
        Mateps_true(:,jj)=Mateps_true(:,jj)+((MatAGE_t*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+MatAGE_t*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    end
    %Last quantile.
    Mateps_true(:,jj)=Mateps_true(:,jj)+(MatAGE_t*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
    
    % Laplace tails
    Mateps_true(:,jj)=Mateps_true(:,jj)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
    
    
    % Earnings
    
    Y(:,jj)=Mateta_true(:,jj)+Mateps_true(:,jj);
    
    
    
    
    % Assets and consumption
    
    
    XX=[];
    for mm1=0:M1
        for mm2=0:M2
            for mm3=0:M3
                for mm4=0:M4
                    XX=[XX hermite(mm1,(A(:,jj)-meanA)/stdA).*hermite(mm2,(Mateta_true(:,jj)-meanY)/stdY).*hermite(mm3,(Y(:,jj)-Mateta_true(:,jj)-meanY)/stdY).*hermite(mm4,(aa-meanAGE)/stdAGE)];
                end
            end
        end
    end
    
    C(:,jj)=XX*Restrue(1:(M1+1)*(M2+1)*(M3+1)*(M4+1),1)+sqrt(Restrue((M1+1)*(M2+1)*(M3+1)*(M4+1)+1))*randn(N,1);
    
    
    
    if jj<=nage-1
        
        
%         if jj==1
%             E_draw=tau3*ones(N,1);
%         else
%             E_draw=randn(N,1);
%         end
%         

        %E_draw=randn(N,1);
        %A(:,jj+1)=[ones(N,1) A(:,jj) Y(:,jj) C(:,jj)]*Res_assets+sqrt(sig_assets)*E_draw;
        
       A(:,jj+1)=log(max((1+.03)*exp(A(:,jj))+exp(Y(:,jj))-exp(C(:,jj)),ones(N,1)));
         
        Mat=zeros(N,(K1+1)*(K2+1));
        for kk1=0:K1
            for kk2=0:K2
                Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true(:,jj)-meanY)/stdY).*hermite(kk2,((aa+2)-meanAGE)/stdAGE);
            end
        end
        
        % age aa_ref+2 receives shock tau0
        if jj==1
            V_draw=tau0*ones(N,1);
        else
            V_draw=unifrnd(0,1,N,1);
        end
        
        %First quantile
        
        Mateta_true(:,jj+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1));
        for jtau=2:Ntau
            Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+...
                ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
                (Vectau(jtau)-Vectau(jtau-1)).*...
                (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
                (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
        end
        %Last quantile.
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+(Mat*Resqtrue(:,Ntau)).*...
            (V_draw>Vectau(Ntau));
        
        % Laplace tails
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
            -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
        
    end
end


figure
plot((aa_ref:2:63)',[mean(Y)'-mean(Y1)'])
figure
plot((aa_ref:2:63)',[mean(C)'-mean(C1)'])
figure
plot((aa_ref:2:63)',[mean(A)'-mean(A1)'])

fig_Y=mean(Y)'-mean(Y1)';
fig_C=mean(C)'-mean(C1)';



% Percentile initial assets
tau3=.90;
V_draw=tau3;

% Proposal, a1


%First quantile
A(:,1)=(Mateta1*Restrue_a1(:,1)).*(V_draw<=Vectau(1));
for jtau=2:Ntau
    A(:,1)=A(:,1)+((Mateta1*(Restrue_a1(:,jtau)-Restrue_a1(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
        (V_draw-Vectau(jtau-1))+Mateta1*Restrue_a1(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
end
%Last quantile.
A(:,1)=A(:,1)+(Mateta1*Restrue_a1(:,Ntau)).*(V_draw>Vectau(Ntau));

% Laplace tails
A(:,1)=A(:,1)+((1/(b1true_a)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
    -(1/bLtrue_a*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));




%%%% FIRST SIMULATION

% percentile of initial shock
tau0=.50;

for jj=1:nage
    
    aa=aa_ref+(jj-1)*2;
    
    MatAGE1=[];
    for kk3=0:K3
        MatAGE1=[MatAGE1 hermite(kk3,(aa-meanAGE)/stdAGE)];
    end
    
    % epsilons
    
    MatAGE_t=[];
    for kk4=0:K4
        MatAGE_t=[MatAGE_t hermite(kk4,(aa-meanAGE)/stdAGE)];
    end
    
    
    
    % Proposal, eta_0
    V_draw=unifrnd(0,1,N,1);
    
    %First quantile
    Mateps_true(:,jj)=(MatAGE_t*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
    for jtau=2:Ntau
        Mateps_true(:,jj)=Mateps_true(:,jj)+((MatAGE_t*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+MatAGE_t*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    end
    %Last quantile.
    Mateps_true(:,jj)=Mateps_true(:,jj)+(MatAGE_t*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
    
    % Laplace tails
    Mateps_true(:,jj)=Mateps_true(:,jj)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
    
    
    % Earnings
    
    Y(:,jj)=Mateta_true(:,jj)+Mateps_true(:,jj);
    
    
    
    
    % Assets and consumption
    
    
    XX=[];
    for mm1=0:M1
        for mm2=0:M2
            for mm3=0:M3
                for mm4=0:M4
                    XX=[XX hermite(mm1,(A(:,jj)-meanA)/stdA).*hermite(mm2,(Mateta_true(:,jj)-meanY)/stdY).*hermite(mm3,(Y(:,jj)-Mateta_true(:,jj)-meanY)/stdY).*hermite(mm4,(aa-meanAGE)/stdAGE)];
                end
            end
        end
    end
    
    C(:,jj)=XX*Restrue(1:(M1+1)*(M2+1)*(M3+1)*(M4+1),1)+sqrt(Restrue((M1+1)*(M2+1)*(M3+1)*(M4+1)+1))*randn(N,1);
    
    
    
    if jj<=nage-1
        
        
        %A(:,jj+1)=[ones(N,1) A(:,jj) Y(:,jj) C(:,jj)]*Res_assets+sqrt(sig_assets)*randn(N,1);
      A(:,jj+1)=log(max((1+.03)*exp(A(:,jj))+exp(Y(:,jj))-exp(C(:,jj)),ones(N,1)));
         
        Mat=zeros(N,(K1+1)*(K2+1));
        for kk1=0:K1
            for kk2=0:K2
                Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true(:,jj)-meanY)/stdY).*hermite(kk2,((aa+2)-meanAGE)/stdAGE);
            end
        end
        
        % age 37 receives shock tau0
        if jj==1
            V_draw=tau0*ones(N,1);
        else
            V_draw=unifrnd(0,1,N,1);
        end
        
        %First quantile
        
        Mateta_true(:,jj+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1));
        for jtau=2:Ntau
            Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+...
                ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
                (Vectau(jtau)-Vectau(jtau-1)).*...
                (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
                (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
        end
        %Last quantile.
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+(Mat*Resqtrue(:,Ntau)).*...
            (V_draw>Vectau(Ntau));
        
        % Laplace tails
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
            -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
        
    end
end

Y1=Y;
C1=C;
Mateta_true1=Mateta_true;
A1=A;

%%%% SECOND SIMULATION

% percentile of initial shock
tau0=.10;

for jj=1:nage
    
    aa=aa_ref+(jj-1)*2;
    
   
    
    MatAGE1=[];
    for kk3=0:K3
        MatAGE1=[MatAGE1 hermite(kk3,(aa-meanAGE)/stdAGE)];
    end
    
    % epsilons
    
    MatAGE_t=[];
    for kk4=0:K4
        MatAGE_t=[MatAGE_t hermite(kk4,(aa-meanAGE)/stdAGE)];
    end
    
    
    
    % Proposal, eta_0
    V_draw=unifrnd(0,1,N,1);
    
    %First quantile
    Mateps_true(:,jj)=(MatAGE_t*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
    for jtau=2:Ntau
        Mateps_true(:,jj)=Mateps_true(:,jj)+((MatAGE_t*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
            (V_draw-Vectau(jtau-1))+MatAGE_t*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
    end
    %Last quantile.
    Mateps_true(:,jj)=Mateps_true(:,jj)+(MatAGE_t*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
    
    % Laplace tails
    Mateps_true(:,jj)=Mateps_true(:,jj)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
        -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
    
    
    
    % Earnings
    
    Y(:,jj)=Mateta_true(:,jj)+Mateps_true(:,jj);
    
    
    
    
    % Assets and consumption
    
    
    XX=[];
    for mm1=0:M1
        for mm2=0:M2
            for mm3=0:M3
                for mm4=0:M4
                    XX=[XX hermite(mm1,(A(:,jj)-meanA)/stdA).*hermite(mm2,(Mateta_true(:,jj)-meanY)/stdY).*hermite(mm3,(Y(:,jj)-Mateta_true(:,jj)-meanY)/stdY).*hermite(mm4,(aa-meanAGE)/stdAGE)];
                end
            end
        end
    end
    
    C(:,jj)=XX*Restrue(1:(M1+1)*(M2+1)*(M3+1)*(M4+1),1)+sqrt(Restrue((M1+1)*(M2+1)*(M3+1)*(M4+1)+1))*randn(N,1);
    
    
    
    if jj<=nage-1
        
        
%         if jj==1
%             E_draw=tau3*ones(N,1);
%         else
%             E_draw=randn(N,1);
%         end
%         

        %E_draw=randn(N,1);
        %A(:,jj+1)=[ones(N,1) A(:,jj) Y(:,jj) C(:,jj)]*Res_assets+sqrt(sig_assets)*E_draw;
        
        A(:,jj+1)=log(max((1+.03)*exp(A(:,jj))+exp(Y(:,jj))-exp(C(:,jj)),ones(N,1)));
         
        Mat=zeros(N,(K1+1)*(K2+1));
        for kk1=0:K1
            for kk2=0:K2
                Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true(:,jj)-meanY)/stdY).*hermite(kk2,((aa+2)-meanAGE)/stdAGE);
            end
        end
        
        % age aa_ref+2 receives shock tau0
        if jj==1
            V_draw=tau0*ones(N,1);
        else
            V_draw=unifrnd(0,1,N,1);
        end
        
        %First quantile
        
        Mateta_true(:,jj+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1));
        for jtau=2:Ntau
            Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+...
                ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
                (Vectau(jtau)-Vectau(jtau-1)).*...
                (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
                (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
        end
        %Last quantile.
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+(Mat*Resqtrue(:,Ntau)).*...
            (V_draw>Vectau(Ntau));
        
        % Laplace tails
        Mateta_true(:,jj+1)=Mateta_true(:,jj+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
            -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
        
    end
end


figure
plot((aa_ref:2:63)',[mean(Y)'-mean(Y1)'])
figure
plot((aa_ref:2:63)',[mean(C)'-mean(C1)'])
figure
plot((aa_ref:2:63)',[mean(A)'-mean(A1)'])

fig_Y1=mean(Y)'-mean(Y1)';
fig_C1=mean(C)'-mean(C1)';

close all

figure 
plot((aa_ref:2:63)',[fig_Y])
figure 
plot((aa_ref:2:63)',[fig_C fig_C1])


