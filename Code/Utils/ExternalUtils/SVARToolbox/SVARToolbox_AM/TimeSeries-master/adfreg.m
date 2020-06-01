function [beta, se, t, tsig, resid, rss, sigma] = adfreg (series, dlags)
%ADFREG (Augmented) Dickey-Fuller regression
%
% [BETA, SE, T, TSIG, RESID, RSS, SIGMA] = ADFREG (SERIES, DLAGS) runs and evaluates a
% Dickey-Fuller (1979) type OLS regression on SERIES with the number of augmented terms
% (differenced lags) given by DLAGS (default is 0):
%                                                                                + resid;
% dseries = beta_0 + beta_1*series(-1) [+ beta_2*dseries(-1) + beta_3*dseries(-2) + ...]
%
% returning: beta = the estimated regression coefficients [beta_0; beta_1; ...]
%              se = the estimated standard errors corresponding to beta
%               t = the t-ratios computed from beta and se
%            tsig = the level at which the corresponding t-ratios are statistically
%                   significantly, using Dickey-Fuller critical values for beta_1 and
%                   standard t-table critical values for beta_2, beta_3, etc.
%                   The significance level of the constant beta_0 is usually of little
%                   interest and since it follows neither a t-distribution nor a Dickey-
%                   Fuller distribution), it is not evaluated here and "NaN" is returned.
%           resid = the (vector) of residuals
%             rss = the residual sum of squares
%           sigma = the estimated standard error of the residuals
%
% (So, reject a unit root if t_1 in the (A)DF regression is statistically significant,
%  e.g. tsig_1 <= 0.05, AND the residuals are not correlated (use, for example, the
%  one of the Durbin or Box-Pierce statistics to check this), because otherwise the test
%  is inefficient.)
%
% REQUIRES the MATLAB Statistics Toolbox and L. Kanzler's DFCRIT m-function.
%
%   The author assumes no responsibility for errors or damage resulting from usage. All
%   rights reserved. Usage of the programme in applications and alterations of the code
%   should be referenced. This script may be redistributed if nothing has been added or
%   removed and nothing is charged. Positive or negative feedback would be appreciated.

%                     Copyright (c) 6 April 1998  by Ludwig Kanzler
%                     Department of Economics, University of Oxford
%                     Postal: Christ Church,  Oxford OX1 1DP,  U.K.
%                     E-mail: ludwig.kanzler@economics.oxford.ac.uk
%                     Homepage:      http://users.ox.ac.uk/~econlrk
%                     $ Revision: 1.01 $      $ Date: 20 May 1998 $

% Compute the total no. of observations in the sample & the series of first differences:
series  = series(:);
obs     = length(series);
dseries = series(2:obs) - series(1:obs-1);

% Create the matrices for the dependent and independent variables:
X     = [series(dlags+1:obs-1)];
if dlags  
    for i = 1 : dlags
        X = [X dseries(dlags+1-i:obs-1-i)];  
    end
end
X     = [ones(obs-dlags-1,1) X];
y     = dseries(dlags+1:obs-1);
% "Run" the OLS REGRESSION, thus obtaining:
beta  = X \ y;                       % vector of est. slope coeff.s (const., lag [dlags])
resid = y - X*beta ;                 % vector of residuals
rss   = resid'*resid;                % residual sum of squares
sgma2 = rss / (obs-2-dlags);         % estimated variance of the residuals, sigma squared
sigma = sqrt(sgma2);                 % ... and sigma
se    = sqrt(diag(sgma2*inv(X'*X))); % vector of est. standard errors of est. coeff.s
t     = beta ./ se;                  % vector of corresponding t-ratios
% Evaluate the SIGNIFICANCE OF ALL T-RATIOS - SEE THE AUTHOR'S SEPARATE DFCRIT.M SCRIPT:
tsig(1:2,1) = [NaN; dfcrit(t(2), obs)];
if dlags > 0
   tsig(3:2+dlags) = min(1 - tcdf(t(3:end), obs-1-dlags), tcdf(t(3:end), obs-1-dlags))*2;
end

% End of function.

% REFERENCE:
%
% Dickey, David & Wayne Fuller (1979), "Distribution of the Estimators for Autoregressive
%    Time Series With a Unit Root", Journal of the American Statistical Association, vol.
%    74, no. 366 (June), pp. 427-431

% End of file.
