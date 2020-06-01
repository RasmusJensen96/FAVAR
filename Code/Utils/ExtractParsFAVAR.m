function [Lag1Pars] = ExtractParsFAVAR(X_st, Y, plag, K, det, slowcode,Dates, WinSize, StartOoS,type)
% Refits the FAVAR and forecasts one-period ahead in a rolling window style
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
%         type     = Boolean factor extaction proc to use, PCA=1 or Kalman = 2 
%------------------------------Outputs-------------------------------------
% Rasmus M. Jensen 2020

Window = StartOoS-WinSize-1:StartOoS-1;
t      = StartOoS; 

X_st_Train = X_st(Window,:); Y_Train = Y(Window,:);
%Dates_Train = Dates(Window);

Date = Dates(t);

[FAVAR, ~] = EstimateReducedFAVAR(X_st_Train, Y_Train, plag, K, det, slowcode,type);
FT   = FAVAR.Ft;

Lag1Pars = FT(1:6,:);
end

