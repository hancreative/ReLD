function [ Shat_mod, Bg ] = ReProCS_final( I, imSize )

Kmin     = 3;
Kmax     = 10;
alpha    = 20;

global rhat;
p = size(I,1);
q = size(I,2);

%training with PCP
mu0         = mean(I,2);
M           = I-repmat(mu0,1, q);
numTrain   = min(100,q);
DataTrain  = M(:,1:numTrain);
lambda = 1/sqrt(max(size(DataTrain)));
[tempL,~,~] = inexact_alm_rpca(DataTrain, lambda);

MTrain      = tempL;
 
[U, Sig, ~] = svd(1/sqrt(numTrain)*MTrain,0);
evals       = diag(Sig).^2;
energy      = sum(evals);
cum_evals   = cumsum(evals);
ind0        = find(cum_evals < 0.95*energy);
rhat        = min(length(ind0),round(numTrain/10));
if rhat<5
 rhat=5;
end

U0          = U(:, 1:rhat); 
lambda_min = evals(rhat);

Shat_mod    = zeros(p,q);
Shat_mod_thresh = zeros(p,q);
Lhat_mod    = zeros(p,q); 
Nhat_mod    = cell(q,1); 
Ut=U0;

Tpred       = [];
Lhat_mod(:,1:numTrain)=MTrain;
Bg(:,1:numTrain) = MTrain+repmat(mu0,1,numTrain);
Shat_mod(:,1:numTrain) = I(:,1:numTrain)-Bg(:,1:numTrain);


Pstar    = Ut;
k        = 0;
K        = [];
addition = 0;
cnew     = [];
t_new    = []; time_upd = []; thresh_diff = []; thresh_base = [];

for t = numTrain+1: size(I,2)
  
clear opts; 
opts.tol   = 1e-3; 
opts.print = 0;
Atf.times  = @(x) Projection(Ut,x); Atf.trans = @(y) Projection(Ut,y);
yt         = Atf.times(M(:,t));


        opts.delta =   norm(Atf.times(Lhat_mod(:, t-1)));
  

    [xp,~]   = yall1(Atf, yt, opts); % xp = argmin ||x||_1 subject to ||yt - At*x||_2 <= opts.delta
    omega(t) = 3*sqrt(M(:,t)'*M(:,t)/p);
    That     = find(abs(xp)>=omega(t));


    Shat_mod_thresh(That,t)    = subLS(Ut,That,yt);
    Shat_mod(:,t) = xp;
    Lhat_mod(:,t) = M(:,t) - Shat_mod(:,t);
    Lhat_mod_thresh(:,t)       = M(:,t) - Shat_mod_thresh(:,t);
    %Fg(That,t)          = I(That,t);
    Nhat_mod{t}         = That;
    Tpred               = That;
    Bg(:,t)             = Lhat_mod(:,t) + mu0;

    
    
%    Projection PCA   
    if addition==0    %&& norm( (Lhat(:,t-alpha+1:t) - Phat*(Phat'*Lhat(:,t-alpha+1:t)))./sqrt(alpha) )>thresh
        addition = 1;
        t_new    = t;
        Pnewhat  = [];
        k        = 0;
    end
        
    if addition==1&&mod(t-t_new+1,alpha)==0
        time_upd = [time_upd,t];           
        D        = Lhat_mod_thresh(:,t-alpha+1:t)-Pstar*(Pstar'*Lhat_mod_thresh(:,t-alpha+1:t)); 
       
        [Pnew_hat, Lambda_new,~] = svd(D./sqrt(alpha),0);
        Lambda_new               = diag(Lambda_new).^2;
        Lambda_new               = Lambda_new(Lambda_new>=lambda_min);
        th                       = round(rhat/3);
        if size(Lambda_new,1)> th
            Lambda_new=Lambda_new(1:th);
        end
           if numel(Lambda_new)==0
               addition  = 0; 
               cnew      = [cnew 0];
           else              
               cnew_hat    = numel(Lambda_new);
               Pnewhat_old = Pnewhat;
               Pnewhat     = Pnew_hat(:,1:cnew_hat); cnew = [cnew cnew_hat];%size(Pnewhat,2)];
               Ut          = [Pstar Pnewhat];   
               
               k=k+1;
               
               if k==1 
                   temp        =(Pnewhat*(Pnewhat'*Lhat_mod_thresh(:,t-alpha+1:t)));
                   thresh_base = [thresh_base norm(temp./sqrt(alpha))];
                   thresh_diff = [thresh_diff norm(temp./sqrt(alpha))];                 
               else
                   temp        =(Pnewhat*(Pnewhat'*Lhat_mod_thresh(:,t-alpha+1:t)));
                   thresh_base = [thresh_base norm(temp./sqrt(alpha))];
                   
                   temp        = (Pnewhat*(Pnewhat'*Lhat_mod_thresh(:,t-alpha+1:t)) - Pnewhat_old*(Pnewhat_old'*Lhat_mod_thresh(:,t-alpha+1:t)));
                   thresh_diff = [thresh_diff norm(temp./sqrt(alpha))];  
               end
               
               flag = 0;
               if k >= Kmin
                   numK = 3;
                   flag = thresh_diff(end)/thresh_base(end-1)<0.01;
                   for ik = 1:numK-1
                       flag = flag && thresh_diff(end-ik)/thresh_base(end-ik-1)<0.01;
                   end
               end
               
               if  k==Kmax|| (k>=Kmin && flag==1)                  
                   addition =0;
                   K        = [K k];
                   Pstar    = Ut;            
               end
           end
    end
    
    
    
    
 
end



end