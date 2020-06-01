function [FE, Fit, Obs, Date] = ExpandingWindowFAVAR(X_st, Y, plag, K, det, slowcode,Dates, StartOoS,  type)
% Refits the FAVAR and forecasts one-period ahead in an expanding window style
% For iterative study use the function in a for-loop
%-------------------------------Inputs-------------------------------------
%         X_st     = Factor data
%         Y z      = VAR data
%         plag     = lags
%         K        = Number of factors 
%         det      = determistics, 0 = none, 1 = const, 2 = lintrend...
%         Dates    = Date-vector
%         slowcode = Factor "speed"
%         StartOoS = Sample split for training and test
%         Var      = Variable to predict, 1:6 in our case (F, Y) %% Removed
%-------------------------------Outputs------------------------------------
%         FE       = One-period ahead forecast error
%         Fit      = Predicted value (t+1)
%         Obs      = Observed value (t+1)
%         Date     = Forecasted date

% Rasmus M. Jensen 2020

t   = StartOoS; % First OoS obs.
tm1 = t-1;      % Last IS obs.

X_st_Train = X_st(1:tm1,:); Y_Train = Y(1:tm1,:);
Dates_Train = Dates(1:tm1);

Date = Dates(t);

[FAVAR, FAVARopt] = EstimateReducedFAVAR(X_st_Train, Y_Train, faclag, plag, K, det, slowcode, type);

Endo = [];
%for i=1:plag
for i=0:plag-1
   Endo = [Endo, FAVAR.ENDO(end-i,:)];
end

FT   = FAVAR.Ft;
Fit  = Endo * FT;

FE  = Y(t,:) - Fit(:,K+1:end); % Residual
Obs = Y(t,:);                  % Observed
end

