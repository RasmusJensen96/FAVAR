function [NumLags, IC_table] = DetermineLagLengthFAVAR(NumFacs,faclag, MaxLag, Y, X_st, slowcode, det,type);
% Inputs X = standardized X
%        Y = VAR Data
%        MaxFac/MaxLag maximum allowed dimensions.
%        faclag = lags in measurement equation
%        slowcode  = Binary indicating slow moving vars in X
%        det = deterministics in the VAR-part: 0 = none, 1 = const, 2 =
%              const + linear trend, 3 = const + quadratic trend.
% Output: Output optimum dimension

% Rasmus M. Jensen 2020


% Number of factors & lags:
K = NumFacs;               % Number of Factors       
% Extract principal components from X


[F0] = DFM_extraction(X_st, faclag, K, type);
% Now rotate the factor space as in Bernanke, Boivin and Eliasz (2005)
% slowindex = find(slowcode==1)';
% xslow = X_st(:,slowindex);
% [Fslow0] = DFM_extraction(xslow, faclag, K, type);
% Fr0 = facrot(F0,Y(:,end),Fslow0);
FY=[F0,Y];  % the extracted factors, plus Y (infl,unemp and interest)

[NumLags, AIC, SCC, HQC,FPE,logL] = VAR_IC(FY,MaxLag,det);

Lags = 1:MaxLag;
IC_table = table(Lags', AIC, SCC, HQC,FPE, logL);
IC_table.Properties.VariableNames{end} = 'Log Likelihood';
IC_table.Properties.VariableNames{1} = 'Lags';
table2latex(IC_table);
end
