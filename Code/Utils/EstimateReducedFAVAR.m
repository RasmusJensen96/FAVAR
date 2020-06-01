function [FAVAR, FAVARopt] = EstimateReducedFAVAR(X, Y, faclag, plag, NumFacs, det, slowcode, type, namesXY,ML, Chol,Include_Y,FC)
%{                         
Estimates and fits an FAVAR with dynamic factor lags:
--------------------------------Inputs-------------------------------------
        X = Factor data (Standardized)
        Y = Observable factor data
        plags   = Number of lags in the FAVAR
        faclags = Number of lags in the DFM
        NumFacs = Number of factors to include
        det = Deterministics in the FAVAR, 0 = none
                                           1 = constant
        slowcode = Binary vector of fast/slow specifications for X
        type     = Boolean switch 1 =  PCA-factors, 2 = Kalman smoothed
                   factors
        NamesXY  = Names of data in the Panel [X', Y']'
        ML       = 1 = DCC-GARCH(1,1) numerical maximum likelihood
                   0 = Analytical Gaussian OLS; eps =dist= N(0,1)
Chol     = Switch determining ordering estimation and factor estimation proc. 
        IncludeY = Switch include y in measurement equation 
                (0): X_t = Lamba(L)*Ft+ eta_t
                (1): X_t = Lamba_F(L)*F_t + Lambda_Y(L)Y_t + eta_t
        FC       = Optional input. Used as a switch with
                   RollingWindowFAVARcomp.m. But should in general not be
                   included or set at 0.
--------------------------------Outputs------------------------------------
        FAVAR    = Fitted FAVAR model with reduced form estimates and DFM
                   as entries
        FAVARopt = Struct containing auxiliary options from the FAVAR
                   estimation for further analysis with the SVAR-toolbox

Rasmus M. Jensen 2020
%}
if ~exist('FC','var')
   FC = 0; 
end
if ~exist('Include_Y','var')
    Include_Y = 0;
end
if ~exist('ML','var') % Checks whether cond. VAR is declared, if not assume IID WN residuals
   ML = 0; 
end

K = NumFacs;               % Number of Factors
%----------------------------Dynamic Factor Model--------------------------

if Include_Y == 0
    [F0] = DFM_extraction(X, faclag, K, type);
else
    [F0] = DFM_extractionWY(X,Y, faclag, K, type);
    F0 = F0(:,1:K);
end
T = size(X,1);
entFFR     = find(strcmp(namesXY,'FFR')) - size(X,2); %% Index to rotate factors around Bernanke ident
if FC == 0
% Rotating factor space as in Bernanke, Boivin and Eliasz (2005)
if Chol == 1
    slowindex = find(slowcode==1)';
    xslow = X(:,slowindex);
    if Include_Y == 0
        [Fslow0] = DFM_extraction(xslow, faclag, K, type);
    else
        [Fslow0] = DFM_extractionWY(xslow,Y, faclag, K, type);
         Fslow0 = Fslow0(:,1:K);
    end
    Fr0 = facrot(F0,Y(:,entFFR),Fslow0);
else
    Fr0 = F0;
end
%-------------------------------FAVAR Model--------------------------------
% Put it all in state-space representation, write obs equ as XY=FY*L+e
XY=[X,Y];    
FY=[Fr0,Y]; 
if Chol == 0
   entCPI     = find(strcmp(namesXY, 'CPI')) - size(X,2); %% Index to CPI
   entIP      = find(strcmp(namesXY, 'IP')) - size(X,2); %% Index to IP
   F1 = 1;
   F2 = 2:NumFacs;
   BQordering = [F1, NumFacs+entIP,NumFacs+entFFR,F2,NumFacs+entCPI];
   FY = FY(:,BQordering);
end
else
    FY = [F0, Y];
    Fr0 = F0;
    XY = [X,Y];
end

[FAVAR, FAVARopt] = VARmodel(FY,plag,det,ML,1); %% Fit FAVAR with ML

XY1 = XY(faclag:end,:);
FYl = FY(faclag:end,:);
for ii=1:faclag-1
     FYl = [FYl, FY(faclag-ii:end-ii,:)];
end
DynamicL = (olssvd(XY1,FYl))';

e = var(XY1 - FYl*DynamicL');

vnames = [];
namei  = 1;
for i = 1:size(XY1,2)
    Xfoo =  XY1(:,i);
    if faclag > 1;
    diffF = size(FYl,2)/size(FY,2) + size(Fr0,2);
    factorindex = 1:diffF:size(FYl,2);
    else
        factorindex = 1;
    end
        for j = 1:size(FY,2)
            Yfoo = FYl(:,factorindex);
            fooreg = nwest(Xfoo,Yfoo,0);
            mR2s(i,j) = fooreg.rsqr;
            if size(factorindex,2) > 1
                indicator = 1;
                fooreg1l = nwest(Xfoo,Yfoo(:,1),0);
                mR21Lags(i,j) = fooreg1l.rsqr;
            end
            factorindex = factorindex + 1;
            if i == size(XY1,2);
                vnames = [vnames; {['Factor ' num2str(namei)]}];
                namei = namei + 1; 
            end
        end
end
vnames(end-size(Y,2)+1:end) = namesXY(end-size(Y,2)+1:end);
if FC == 0
    if Chol == 0
        vnames = vnames(BQordering);
    end
end
tablemR2 = array2table(mR2s); tablemR2.Properties.RowNames = namesXY; tablemR2.Properties.VariableNames = vnames;
if size(factorindex,2) > 1
tablemR21lag = array2table(mR21Lags); tablemR21lag.Properties.RowNames = namesXY; tablemR21lag.Properties.VariableNames = vnames;
end


SIGMA = diag(diag(e'*e./T));

L = (olssvd(XY,FY))';
R2s = 1 - e; 

R2s = table(R2s'); R2s.Properties.RowNames = namesXY; R2s.Properties.VariableNames = "$R^2$";

V1   = sum(diag(cov(X)))-size(Y,2);
    VDFM = sum(diag(cov(FYl*DynamicL')))-size(Y,2);
    V    = VDFM/V1;
    VDFMs = sum(diag(cov(FY*L')))-size(Y,2);
    Vs    = VDFMs/V1;
FAVAR.DFM.Vstatic = Vs;

if ML == 1
    FAVAR.GARCHresiduals = true;
else
    FAVAR.GARCHresiduals = false;
end
FAVAR.NumberFactors = NumFacs;
FAVAR.DFM.L = L;
FAVAR.DFM.DynamicL = DynamicL;
FAVAR.DFM.mR2s = tablemR2;
if exist('indicator','var')
    FAVAR.DFM.mR21lags = tablemR21lag;
end
FAVAR.DFM.R2s = R2s;
FAVAR.DFM.IdiosyncErr = e;
FAVAR.DFM.ENDOx = XY;
FAVAR.DFM.ENDOy = FY;
FAVAR.DFM.SIGMA = diag(diag(SIGMA));
FAVAR.DFM.Am1  = chol(SIGMA);
FAVAR.DFM.Am1  = FAVAR.DFM.Am1';
    FAVAR.DFM.V = V;
FAVAR.DFM.FactorLags = faclag;

foo = round(var(Y),1);
if foo(1)==1 && foo(2)==1 && foo(3)==1
    FAVARopt.standardized = true;
else
    FAVARopt.standardized = false;
end

end