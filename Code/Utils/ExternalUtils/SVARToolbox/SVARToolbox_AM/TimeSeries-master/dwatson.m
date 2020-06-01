function [dw, dwsigup, dwsiglow] = dwatson (series)
%DWATSON Durbin-Watson d-statistic and significance level for the null hypothesis: DW = 2
%
% [DW, DWSIGUP, DWSIGLOW] = DWATSON(SERIES) computes
%
%    - the Durbin-Watson (1950, 1951) statistic d of serial correlation, and
%
%    - the significance level, if any, at which the null hypothesis DW = 2 is rejected
%      against either of the one-sided alternatives (but not both!), using both upper-
%      bound and lower-bound critical values. Possible values are 0.01, 0.05, 1, NaN.
%
%      This lookup is valid only for regressions with one explanatory variable & constant.
%
% The author assumes no responsibility for errors or damage resulting from usage. All
% rights reserved. Usage of the programme in applications and alterations of the code
% should be referenced. This script may be redistributed if nothing has been added or
% removed and nothing is charged. Positive or negative feedback would be appreciated.

%                     Copyright (c) 17 March 1998 by Ludwig Kanzler
%                     Department of Economics, University of Oxford
%                     E-mail:       ludwig.kanzler-alumni@lse.ac.uk
%                     Homepage:   http://www2.gol.com/users/kanzler
%                     $ Revision: 1.3 $        $ Date: 5 May 2005 $

% Compute the Durbin-Watson statistic, which is defined as
% sum((series(2:end)-series(1:end-1)).^2)/sum((series(2:end)-mean(series(2:end))).^2),
% but which can be computed much quicker as follows:
series  = series(:);
dseries = series(2:end)-series(1:end-1);
dw      = diag((dseries'*dseries)./(series'*series))';

if nargout > 1
   % DWCRIT lists the critical values for the d-statistic for one explanatory variable
   % (excluding the constant) as published in Durbin & Watson (1951) and Savin & White
   % (1977), reproduced in Judge et al. (1988). Add entries if you need/know more values.
   %
   %                                  sample
   %                                   size
   dwcrit =  [  NaN   NaN  0.610 1.400   6
                NaN   NaN  0.700 1.356   7
                NaN   NaN  0.763 1.332   8
                NaN   NaN  0.824 1.320   9
                NaN   NaN  0.879 1.320  10
                NaN   NaN  0.927 1.324  11
                NaN   NaN  0.971 1.331  12
                NaN   NaN  1.010 1.340  13
                NaN   NaN  1.045 1.350  14
               0.81  1.07  1.077 1.361  15
               0.84  1.09  1.106 1.371  16
	       0.87  1.10  1.133 1.381  17
               0.90  1.12  1.158 1.391  18
               0.93  1.13  1.180 1.401  19
               0.95  1.15  1.201 1.411  20
               0.97  1.16  1.221 1.420  21
               1.00  1.17  1.239 1.429  22
               1.02  1.19  1.257 1.437  23
               1.04  1.20  1.273 1.446  24
               1.05  1.21  1.288 1.454  25
               1.07  1.22  1.302 1.461  26
               1.09  1.23  1.316 1.469  27
               1.10  1.24  1.328 1.476  28
               1.12  1.25  1.341 1.483  29
               1.13  1.26  1.352 1.489  30
               1.15  1.27  1.363 1.496  31
               1.16  1.28  1.373 1.502  32
               1.17  1.29  1.383 1.518  33
               1.18  1.30  1.393 1.514  34
               1.19  1.31  1.402 1.519  35
               1.21  1.32  1.411 1.525  36
               1.22  1.32  1.419 1.530  37
               1.23  1.33  1.427 1.535  38
               1.24  1.34  1.435 1.540  39
               1.25  1.34  1.442 1.544  40
               1.29  1.38  1.475 1.566  45
               1.32  1.40  1.503 1.585  50
               1.36  1.43  1.528 1.601  55
               1.38  1.45  1.549 1.616  60
               1.41  1.47  1.567 1.629  65
               1.43  1.49  1.583 1.641  70
               1.45  1.50  1.598 1.652  75
               1.47  1.52  1.611 1.662  80
               1.48  1.53  1.624 1.671  85
               1.50  1.54  1.635 1.679  90
               1.51  1.55  1.645 1.687  95
               1.52  1.56  1.654 1.694 100
                NaN   NaN  1.720 1.746 150
                NaN   NaN  1.758 1.778 200

               0.01  0.01   0.05  0.05 NaN ]; % significance level for
                                              % lower upper lower upper bounds
   [r,c] = size(dwcrit);

   % Find the entry matching sample size and dw statistic:
   matchsize = find(length(series) >= dwcrit(1:r-1, c));
   matchlow  = find(min(dw, 4-dw)  <= dwcrit(matchsize(end),1:2:c-1));
   matchup   = find(min(dw, 4-dw)  <= dwcrit(matchsize(end),2:2:c-1));

   % Assign the corresponding level of significance, if any:
   if length(series) > 1.25*dwcrit(r-1, c)
      dwsiglow = NaN;
      dwsigup  = NaN;
   else   
      if isempty(matchlow)
         dwsiglow = 1;
      else
         dwsiglow = dwcrit(r, fix(1.5*matchlow(1)));
      end
      if isempty(matchup)
         dwsigup = 1;
      else
         dwsigup   = dwcrit(r, 2*matchup(1));
      end

   end
end

% End of function.


% NOTE that there is a case for reporting the significance at the UPPER bound only:
%      the increase in the Type I error arising from ignoring the inconclusive area is
%      likely to be small when (see Harvey, 1990, p. 203)
%      - there is only one dependent variable (so k = 2),
%      - there are a large number of observations, and
%      - the data is in levels).

% THANKS to Daniel Reuman for pointing out a bug in the algorithm for assigning the
%       correct significance level, which has been fixed with this version 1.3.

% REFERENCES:
%
% Durbin, James & G.S. Watson (1950), "Testing for Serial Correlation in Least Squares
%    Regression I", Biometrika, vol. 37, pp. 409-428 [corrections in Durbin & Watson
%    (1951), pp. 177-178]
%
% Durbin, James & G.S. Watson (1951), "Testing for Serial Correlation in Least Squares
%    Regression II", Biometrika, vol. 38, pp. 159-178
%
% Harvey, Andrew (1990), "The Econometric Analysis of Time Series", 2nd edition,
%    Press, Cambridge, Massachusetts, pp. 221-223
%
% Judge, George, Carter Hill, William Griffiths, Helmut LÅEkepohl & Tsoung-Chao Lee
%    (1988), "Introduction to the Theory and Practice of Econometrics", 2nd edition, John
%    Wiley & Sons, New York, pp. 991-995, Table 5, k=2
%
% Savin, Eugene & Kenneth White (1977), "The Durbin-Watson Test for Serial Correlation
%    with Extreme Sample Sizes or Many Regressors", Econometrica, vol. 45, no. 8
%    (November), pp. 1989-1996

% End of file.
