function [Factor, Kprop, A] = DFM_extractionWY(X, Y, plag, K, type)
% Function returning filtered factors from a DFM either with PCA or
% filtered and smoothed with Kalman.
%-------------------------------Inputs-------------------------------------
% X:    Dataset from which the factors should be extracted (t x N)
% Y:    Observable factors
% plag: Number of lags in the DFM; order of Psi(L)
% K:    Number of factors to be extracted
% type: Boolean switch, 1 = PCA, 2 = Kalman Smoothed
%------------------------------Outputs-------------------------------------
% Factors: Extracted timeseries of factors (t x K)
% Rasmus M. Jensen 2020

[t,N] = size(X);

X = (X-mean(X))./std(X);              % Standardize data

%[~,fooval] = eigs(X'*X,size(X,2));       % 
%[eigenvector,Val] = eigs(X'*X,K);        % Extract eigenvectors (Descending order)
[~,fooval] = eigs(cov(X),size(X,2));       % 
[eigenvector,Val] = eigs(cov(X),K);        % Extract eigenvectors (Descending order)

eigenvector = eigenvector*sqrt(N);

TotVar = sum(diag(fooval)); % Total Variance; sum of eigenvalues 

F = X * eigenvector/N;

if type == 1 
    Factor = F;
    Kprop = sum(diag(Val))/sum(diag(fooval));
    A = [];
    return
end
CC = X*eigenvector*eigenvector'/N;      % Common Comp

r = size(Y,2);
F = [F,Y];

A_temp = zeros(K,(K+r)*plag)';
I = eye((K+r)*plag,K*plag);
A = [A_temp';I(1:end-K,1:end)];
Q = eye(plag*K,plag*K);

[F_endo, F_lag] = VARmakexy(F,plag-1,0);  % Generate Reg matrices
L = inv(F_lag'*F_lag) * F_lag'*F_endo;    % VAR-system  
A(1:K,1:K*(plag-1)) = L';
e = F_endo - F_lag*A(1:K,1:K*(plag-1))';  % Factor error
H = e'*e;                                 % Error Comp.
Q(1:K,1:K) = H;                           % Comp.
idiComp = X-CC;                           % Idiosynchratic Comp.
R = diag(diag(idiComp'*idiComp));         % R diagonal
%R = diag(var(idiComp)');
[~,Z] = VARmakexy(F,plag,0);
Xhat_ini = Z(1,:)';                       % Initial state mean
%Xhat_ini = zeros(K*plag, 1);
%PHat_ini = eye(K*plag);
PHat_ini = reshape(pinv(eye((K*plag)^2) - kron(A,A))*vec(Q), K*plag, K*plag);           % Initial state covariance
           % Initial state covariance
C = [eigenvector, zeros(N,K*(plag-1))];
%------------------------------Kalman Calls--------------------------------
[x_filt,x_filtP,Ptt,Pttm]=KalmanFilterDFM(Xhat_ini,PHat_ini,X,A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=KalmanSmootherDFM(A,x_filt,x_filtP,Ptt,Pttm,C,R);
Factor =  xsmooth(1:K,:)';

%------------------------Total Variance Explained--------------------------
[~, foo1val]= eigs(xsmooth*xsmooth');
TotVar1 = sum(diag(foo1val));

Kprop = TotVar1 / TotVar;

% Kprop = sum(diag(Vsmooth(1:K*plag,1:K*plag,1)))/TotVar;

end

function [x_filt,x_Pred,P_filt,P_Pred]=KalmanFilterDFM(initx,initP,x,A,C,R,Q)
%----------------------------------Input-----------------------------------
% x(:,t) = the observation at time t
% A      = the system matrix
% C      = EigVecs init.
% Q      = Filtered factor error ovariance
% R      = Idiosync covariance
% initx  = a priori state mean t=1
% initV  = a priori state cov, t=1
%----------------------------------Output----------------------------------
% x_Pred = E[X(:,t) | y(:,1:t-1)]
% P_Pred = Sig[X(:,t) | y(:,1:t-1)]
% x_filt = E[X(:,t) | y(:,1:t)]
% P_filt = Sig[X(:,t) | y(:,1:t)]

[T,N] = size(x);
r     = size(A,1);
y=x';
x_Pred=[initx zeros(r,T)];
x_filt=zeros(r,T);

P_Pred=zeros(r,r,T);
P_Pred(:,:,1)=initP;
P_filt=zeros(r,r,T);

for j=1:T
    L=inv(C*P_Pred(:,:,j)*C'+R);
    %L = C * (P_Pred(:,:,j) \ C') + R;
    x_filt(:,j)=x_Pred(:,j)+P_Pred(:,:,j)*C' * L * (y(:,j)-C*x_Pred(:,j));
    P_filt(:,:,j)=P_Pred(:,:,j)-P_Pred(:,:,j)*C'*L*C*P_Pred(:,:,j); 
    x_Pred(:,j+1)=A' * x_filt(:,j);
    P_Pred(:,:,j+1)=A'*P_filt(:,:,j)*A+Q;
end
end
function [xitT,PtT,PtTm]=KalmanSmootherDFM(A,x_filt,x_pred,p_filt,P_pred,C,R)
%----------------------------------Input-----------------------------------
% y(:,t) - observation
% A - the system matrix (VAR-estimates)
% xittm = E[X(:,t) | y(:,1:t-1)]
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]
% C - the observation matrix 
% R - the observation covariance
%----------------------------------Output----------------------------------
% xitT = E[X(:,t) | y(:,1:T)]
% PtT = Cov[X(:,t) | y(:,1:T)]
% PtTm = Cov[X(:,t+1),X(:,t) | y(:,1:T)]

[T]=size(x_filt,2);
r=size(A,1);
P_pred=P_pred(:,:,1:end-1);
x_pred=x_pred(:,1:end-1);
J=zeros(r,r,T);
for i=1:T-1
    J(:,:,i)=p_filt(:,:,i)*A'*inv(P_pred(:,:,i+1));
end
for i=1:T
    L(:,:,i)=inv(C*P_pred(:,:,i)*C'+R);
    K(:,:,i)=P_pred(:,:,i)*C'*L(:,:,i);
end
xitT=[zeros(r,T-1)  x_filt(:,T)];
PtT=zeros(r,r,T);
PtTm=zeros(r,r,T);
PtT(:,:,T)=p_filt(:,:,T);
PtTm(:,:,T)=(eye(r)-K(:,:,T)*C)*A*p_filt(:,:,T-1);
for j =1:T-1
    
    xitT(:,T-j)= x_filt(:,T-j)+J(:,:,T-j)*(xitT(:,T+1-j)-x_pred(:,T+1-j));
    PtT(:,:,T-j)=p_filt(:,:,T-j)+J(:,:,T-j)*(PtT(:,:,T+1-j)-P_pred(:,:,T+1-j))*J(:,:,T-j)';
    
    
end
for j =1:T-2
    PtTm(:,:,T-j)=p_filt(:,:,T-j)*J(:,:,T-j-1)'+J(:,:,T-j)*(PtTm(:,:,T-j+1)-A*p_filt(:,:,T-j))*J(:,:,T-j-1)';
end
end