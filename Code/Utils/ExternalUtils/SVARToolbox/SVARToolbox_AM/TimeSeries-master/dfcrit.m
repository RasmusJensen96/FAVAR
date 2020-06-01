function [sig, crit] = dfcrit (tratio, ssize, variant)
%DFCRIT Critical Dickey-Fuller values and level of significance, based on MacKinnon (1991)
%
% [SIG, CRIT] = DFCRIT (TRATIO, SSIZE, VARIANT) computes the critical values CRIT
%    of the Dickey-Fuller distribution for given sample size SSIZE and returns the
%    level SIG, if any, at which t-value TRATIO is significant.
%
%    Critical values are returned as a row vector for the 1%, 5% and 10% significance
%    levels of a one-sided test. SIG is the highest level at which TRATIO is significant;
%    1 indicates that TRATIO was not significant at the 10% level.
%
%    VARIANT should be set to 1, 2 or 3 in accordance to whether the Dickey-Fuller
%    regression contains
%                        (1) no constant and no trend,
%                        (2) a constant but no trend,
%                        (3) a constant and a trend coefficient.
%
%    The relevant critical values CRIT are computed from a response surface developed by
%    MacKinnon (1991).
%
% The author assumes no responsibility for errors or damage resulting from usage. All
% rights reserved. Usage of the programme in applications and alterations of the code
% should be referenced. This script may be redistributed if nothing has been added or
% removed and nothing is charged. Positive or negative feedback would be appreciated.

%                     Copyright (c) 23 April 1998 by Ludwig Kanzler
%                     Department of Economics, University of Oxford
%                     Postal: Christ Church,  Oxford OX1 1DP,  U.K.
%                     E-mail: ludwig.kanzler@economics.oxford.ac.uk
%                     Homepage:      http://users.ox.ac.uk/~econlrk
%                     $ Revision: 1.0 $     $ Date: 23 April 1998 $

if nargin < 3
   variant = 2;
   if nargin < 2
      error('Insufficient number of arguments! Execution aborted.')
   end
elseif ~sum(variant == [1 2 3])
   error('Inadmissable value for VARIANT! Execution aborted.')
end

binf = [  -2.5658   -3.4335   -3.9638
          -1.9393   -2.8621   -3.4126
          -1.6156   -2.5671   -3.1279 ];
   
b1   = [  -1.960    -5.999    -8.353
          -0.398    -2.738    -4.039
          -0.181    -1.438    -2.418  ];
   
b2   = [ -10.04    -29.25    -47.44
           0.0      -8.36    -17.83
           0.0      -4.48     -7.58   ];

crit = binf(:, variant) + b1(:, variant)./ssize + b2(:, variant)./ssize^2;

switch sum(tratio <= crit)
   case 0
      sig = 1;
   case 1
      sig = 0.10;
   case 2
      sig = 0.05;
      if variant == 2
         if dftable (tratio, ssize) == 0.025;
            sig = 0.025;
         end
      end
   case 3
      sig = 0.01;
end

% End of main function.

function dfsig = dftable (tratio, ssize)
%DFTABLE Level of significance for Dickey-Fuller regression with constant, but no trend
%
%                    Copyright (c) 13 March 1998  by Ludwig Kanzler
%                    Department of Economics,  University of Oxford
%                    Postal: Christ Church, Oxford OX1 1DP, England
%                    E-mail:  ludwig.kanzler@economics.oxford.ac.uk
%                    $ Revision: 1.03 $ $ Date: 15 September 1998 $

% DFCRIT lists the critical values for the Phillips-Perron Z(t) Test and for the Dickey-
% Fuller Test Based on Estimated OLS t Statistic, as tabulated in Fuller (1976) on p. 373,
% and reproduced, among others, in Hamilton (1994), p. 763, Table B.6, Case 2, and in
% Harvey (1990), p. 368, Table D, 2nd matrix.
%
% DFCRIT could probably be improved by using the extended tabulations of Guilkey & Schmidt
% (1989). Tables for the cases of no constant and of constant and trend could also be
% added here.
%
%                                   sample
%                                    size
dfcrit = [ -3.75 -3.33 -3.00 -2.63     25
           -3.58 -3.22 -2.93 -2.60     50
           -3.51 -3.17 -2.89 -2.58    100
           -3.46 -3.14 -2.88 -2.57    250
           -3.44 -3.13 -2.87 -2.57    500
           -3.43 -3.12 -2.86 -2.57    inf

            0.01  0.025 0.05  0.10    NaN ];
%          prob(tratio<=table entry)

% Find the entry matching sample size and t-ratio:
matchsize  = find(ssize  >= dfcrit(1:6,5));
matchlevel = find(tratio <= dfcrit(matchsize(end),1:4));

% Assign the corresponding level of significance, if any:
if matchlevel
   dfsig   = dfcrit(7, matchlevel(1));
else
   dfsig   = 1;
end

% End of sub-function.


% REFERENCES:
%
% Fuller, Wayne (1976), "Introduction to Statistical Time Series", John Wiley & Sons, New
%    York
%
% Guilkey, David & Peter Schmidt (1989), "Extended Tabulations for Dickey-Fuller Tests",
%    Economics Letters, vol. 31, no. 4, pp. 355-357
%
% Hamilton, James (1994), "Time Series Analysis", Princeton University Press, Princeton,
%    New Jersey
%
% Harvey, Andrew (1990), "The Econometric Analysis of Time Series", 2nd edition, MIT
%    Press, Cambridge, Massachusetts
%
% MacKinnon, James (1991), "Critical Values for Cointegration Tests", Chapter 13 in
%    Robert Engle & Clive Granger, eds., "Long-run Economic Relationships: Readings
%    in Cointegration", Oxford University Press, Oxford, pp. 267-276

% End of file.
