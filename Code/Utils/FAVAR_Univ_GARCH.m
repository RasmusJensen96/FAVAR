function [Meanmodel, volatilitymodel] = FAVAR_Univ_GARCH(Y,X)
% Function conducting a univariate regression with GARCH(1,1)-residuals

if size(X,1)<size(X,2);
   disp(['Dimensions seems wrong, defined X=t(X)']);
   X = X'; 
end
[T, N] = size(X);

[pars, LLH] = FitLinearModelGARCH(Y,X);

%-------------------------------Mean model---------------------------------
% Regression betas y = xb' +  eps
Meanmodel.betas = pars(1:end-3);                             % Coefficients
Meanmodel.yhat  = X * Meanmodel.betas';                     % Filtered mean
resid   = Y-Meanmodel.yhat;                                     % Residuals
Meanmodel.resid = resid';
Meanmodel.sige = Meanmodel.resid' * Meanmodel.resid/(T-N);% Standard errors
Meanmodel.LLH = LLH;                                            % Total LLH
%-----------------------------Variance model-------------------------------
% Variance parameters s = omega + alpha * eps^2 + beta * sigma^2
omega = pars(end-2); alpha = pars(end-1); beta = pars(end);
[LLH, Sigma2] = GARCHFilter(Meanmodel.resid',omega, alpha, beta);  % Filter
volatilitymodel.filteredvolatility = sqrt(Sigma2);    % Filtered volatility
volatilitymodel.standardizedresidual = Meanmodel.resid'./volatilitymodel.filteredvolatility;
volatilitymodel.LLH = LLH;                        % marginal log-likelihood
volatilitymodel.pars = [omega, alpha, beta];      % parameters
end

