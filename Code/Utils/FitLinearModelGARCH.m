function [pars, LLH] = FitLinearModelGARCH(Y,X)

sizeP = size(X, 2);

initp=X\Y; % Starting values of the regression coefficients
res = Y-X*initp;

Theta = [initp',0,0.1,0.8]; % repmat(0.1,1,3)];
%alpha = 0.1;
%beta = 0.8;
%omega = 1/(1-alpha-beta);
Theta(end-2) = var(res)*(1-Theta(end-1)-Theta(end));%var(res)/(1-Theta(end-1)-Theta(end));

options = optimset('Display','off','TolFun',1e-5,'TolX',1e-5,'MaxFunEvals',40000,'MaxIter',40000);

%[pars, LLH] = fminsearch(@(Theta)NLL_ARDLGARCH(Theta),Theta,options);
[pars, LLH] = fminunc(@(Theta)NLL_ARDLGARCH(Theta),Theta,options);


function LLH = NLL_ARDLGARCH(Theta)
% LLH of a linear regression with gaussian errors

RegPars = Theta(1:end-3);
dOmega = Theta(end-2);
dAlpha = Theta(end-1);
dBeta = Theta(end);

foo = dAlpha + dBeta;
 if foo > 1 || foo<0 || dAlpha>1 || dBeta>1 ||dBeta<0 || dAlpha<0 
       LLH = 1e9; 
       return
 end

u = Y-X * RegPars';
sigma2 = zeros(size(u,1),1);
sigma2(1) = dOmega/(1-dAlpha-dBeta);
LLH = log(normpdf(u(1), 0, sqrt(sigma2(1))));
for i = 2:size(u,1)
   sigma2(i) =  dOmega + dAlpha * u(i-1)^2 + dBeta * sigma2(i-1);
   LLH = LLH + log(normpdf(u(i), 0, sqrt(sigma2(i))));
end
LLH = -LLH;
%LLH = -sum(log(normpdf(u,0, sqrt(sigma2))));

end
end