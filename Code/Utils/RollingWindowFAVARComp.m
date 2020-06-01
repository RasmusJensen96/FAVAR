function [FE, Fit, Obs,Dates_Test] = RollingWindowFAVARComp(X_st, Y,faclag, plag, K, det, slowcode,Dates, WinSize, StartOoS,type, h, singleV,namesXY)
% Refits the FAVAR and forecasts one-period ahead in a rolling window style
% For iterative study use the function in a for-loop
%-------------------------------Inputs-------------------------------------
%         X_st     = Factor data
%         Y z      = VAR data
%         faclag   = Lags in the measurement equation DFM
%         plag     = lags
%         K        = Number of factors 
%         det      = determistics, 0 = none, 1 = const
%         Dates    = Date-vector
%         slowcode = Factor "speed" - % No longer needed
%         WinSize  = Size of training window
%         StartOoS = Sample split for training and test
%         type     = Boolean factor extaction proc to use, PCA=1 or Kalman = 2 
%         h        = Forecast horizon
%         singleV  = Boolean determines whether the function returns a
%                    vector or a single double value t periods ahead.
%         namesXY  = names of the time-series in XY
%------------------------------Outputs-------------------------------------
%         FE       = One-period ahead forecast error y(t+h|t+h)-y(t+h|t)
%         Fit      = Predicted value (t+h|t)
%         Obs      = Observed value (t+h|t+h)
%         Date     = Forecasted date

% Rasmus M. Jensen 2020

Window = StartOoS-WinSize-1:StartOoS-1;
t      = StartOoS; 

X_st_Train = X_st(Window,:); 
Y_Train = Y(Window,:); 
Y_Test  = Y(t-1+h,:);
Dates_Test = Dates(t:t+h-1);

[FAVAR, FAVARopt] = EstimateReducedFAVAR(X_st_Train, Y_Train,faclag, plag, K, det, slowcode,type,namesXY,0,0,0,1);
n = FAVAR.nvar;                       
Endo = [];
%for i=1:plag
for i=0:plag-1
   Endo = [Endo, FAVAR.ENDO(end-i,:)];
end

if det == 1
    %Const = repmat(FAVAR.Ft(1,:), 1,plag);
    Const = [FAVAR.Ft(1,:), zeros(1, (plag-1)*size(FAVAR.ENDO,2))];
else
    Const = 0;
end

Endo = Endo';

FT   = FAVAR.Fcomp;

if singleV == true
    Fcomp_eye = eye(size(FT,1));
    for jj = 1:h
        Fcomp_eye = FT * Fcomp_eye; 
        if det == 1
        C = eye(size(Const, 2)); % Fixing Constant term; Lutkepohl p.34
            for tt = 1:jj-1
               C = C + FT^(tt);
            end
                C  = C*Const';
        else
            C = 0;
        end
        Fit(jj,:)  = C + (Fcomp_eye * Endo)';
        FE( jj,:)  = Y(t+jj-1,:) - Fit(jj,K+1:K+3); % Forecast error
    end
    Dates_Test = [Dates(t-1); Dates(t:t+h-1)];
    Fit = [Y(t-1,:); Fit(:,K+1:K+3)];
    Obs = Y(t-1:t+h-1,:);% Observed
else
iter = 0;
for ii = h
    iter = iter+1;
    if det == 1;
        C = eye(size(Const, 2)); % Fixing Constant term; Lutkepohl p.34
        for jj = 1:ii-1
            C = C + FT^(jj);
        end
        C = C * Const';
    else
        C = 0;
    end
    Fitfoo = C + FT^(ii) * Endo;
    Fitfoo =  Fitfoo(K+1:K+3);
    Fit(:,iter) = Fitfoo; % Extract other observables
    Y_Test = Y(Window(end)+ii,:);
    FE(:,iter)  = Y_Test - Fitfoo';
    Dates_Test(:,iter) = Dates(t+ii-1);
    Obs(:,iter) = Y_Test;              % Observed
    end
    %for i == h
    %Fit = h * constantVector + FT^h * Endo;
    %Fit = Const + (FT^h * Endo)';
    %Fit = Fit(K+1:K+3); % Extract other observables
    %FE  = Y_Test - Fit;
    %Dates_Test = Dates_Test(end);
    %Obs = Y_Test;              % Observed
end

