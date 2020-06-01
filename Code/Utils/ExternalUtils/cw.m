function [cw_stat,cw_pval] = cw(y1,y2,y)
% Purpose: 
% This function implements the Clark-West test for forecasting performance. 
% The null is that MSFE_1 <= MSFE_2, and the alternative is that 
% MSFE_1 > MSFE_2. 
% 
% Inputs:   
% y1: The first forecast.
% y2: The second forecast.
% y:  The original series.
% 
% Outputs:
% cw_stat: The test statistic of the test
% cw_pval: The corresponding p-value of the test statistic.
% 
% Written by:
% Rasmus Fisker Bang, Aarhus University.

% Edit, Rasmus Jensen 2020
if cols(y1) > rows(y1)
    y1 = y1';
end
if cols(y2) > rows(y2)
    y2 = y2';
end
if cols(y) > rows(y)
    y = y';
end
% Defining terms for the CW test statistic
e1     = (y-y1).^2;
e2     = (y-y2).^2;
e3     = (y1-y2).^2;
cwtemp = e1 - e2 + e3;
% Making the regression to get residuals
cwfit_ols = OLSmodel(cwtemp,ones(length(cwtemp),0));

% Checking autocorrelation in the residuals, in order to include the right
% amount of lags in the Newey-West estimator.
[acorr,~,bound] = autocorr(cwfit_ols.resid);
no_lags         = max(min(find((acorr < 3*bound(1) & acorr > 3*bound(2)) == 1))-2,1);
% Making the new regression with Newey-West standard errors.
cwfit_nwest = nwest(cwtemp,ones(length(cwtemp),1),no_lags);
cw_stat     = cwfit_nwest.tstat;
cw_pval     = 1-normcdf(cw_stat,0,1);
end
