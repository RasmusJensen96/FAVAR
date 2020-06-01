function [FAVAR, FAVARopt] = VARmodel(ENDO,nlag,const,ML,DCC,EXOG,nlag_ex)
%{
 Perform Factor augmented vector autogressive (FAVAR) estimation with OLS  
    or Numerical ML of the DCC-GARCH(1,1)-linear regression
-------------------------------------Inputs--------------------------------
	- ENDO: an (nobs x nvar) matrix of y-vectors
	- nlag: lag length
	- const: 0 no constant; 1 constant; 2 constant and trend; 3 constant, 
       trend, and trend^2 [dflt = 0] (Not tested properly for the FAVAR)
    - ML: 1 = Numerical Maximum Likelihood (DCC-GARCH(1,1)-FAVAR), 
          0 = Analytical Gaussian ML N(0,1)-residuals
	- EXOG: optional matrix of variables (nobs x nvar_ex) (Not tested properly in the FAVAR)
       - nlag_ex: number of lags for exogeonus variables [dflt = 0]
------------------------------------Outputs--------------------------------
   - FAVAR: structure including FAVAR estimation results
   - FAVARopt: structure including FAVAR options
---------------------------------------------------------------------------

SVAR toolbox by:
 Ambrogio Cesa Bianchi, March 2015
 ambrogio.cesabianchi@gmail.com
 Note: this code is a modified version of of the vare.m function of James 
 P. LeSage

Modified to accomodate the FAVAR by:
Rasmus M. Jensen 2020
%}
%% Check inputs
%--------------------------------------------------------------------------
[nobs, nvar] = size(ENDO);
if exist('ML','var') == 0 % Checks whether cond. VAR is declared, if not assume IID WN residuals
   ML = 0; 
end
if exist('DCC','var') == 0
   DCC = 0;
end
% Create VARopt and update it
FAVARopt = VARoption;
FAVAR.ENDO = ENDO;
FAVAR.nlag = nlag;

% Check if ther are constant, trend, both, or none
if ~exist('const','var')
    const = 1;
end
FAVAR.const = const;

% Check if there are exogenous variables 
if exist('EXOG','var')
    [nobs2, nvar_ex] = size(EXOG);
    % Check that ENDO and EXOG are conformable
    if (nobs2 ~= nobs)
        error('var: nobs in EXOG-matrix not the same as y-matrix');
    end
    clear nobs2
    % Check if there is lag order of EXOG, otherwise set it to 0
    if ~exist('nlag_ex','var')
        nlag_ex = 0;
    end
    FAVAR.EXOG = EXOG;
else
    nvar_ex = 0;
    nlag_ex = 0;
    FAVAR.EXOG = [];
end


%% Save some parameters and create data matrices
%--------------------------------------------------------------------------
    nobse         = nobs - max(nlag,nlag_ex);
    FAVAR.nobs      = nobse;
    FAVAR.nvar      = nvar;
    FAVAR.nvar_ex   = nvar_ex;    
    FAVAR.nlag      = nlag;
    FAVAR.nlag_ex   = nlag_ex;
    ncoeff        = nvar*nlag; 
    FAVAR.ncoeff    = ncoeff;
    ncoeff_ex     = nvar_ex*(nlag_ex+1);
    ntotcoeff     = ncoeff + ncoeff_ex + const;
    FAVAR.ntotcoeff = ntotcoeff;
    FAVAR.const     = const;

% Create independent vector and lagged dependent matrix
[Y, X] = VARmakexy(ENDO,nlag,const);

% Create (lagged) exogeanous matrix
if nvar_ex>0
    X_EX  = VARmakelags(EXOG,nlag_ex);
    if nlag == nlag_ex
        X = [X X_EX];
    elseif nlag > nlag_ex
        diff = nlag - nlag_ex;
        X_EX = X_EX(diff+1:end,:);
        X = [X X_EX];
    elseif nlag < nlag_ex
        diff = nlag_ex - nlag;
        Y = Y(diff+1:end,:);
        X = [X(diff+1:end,:) X_EX];
    end
end


%% ML estimation equation by equation
%--------------------------------------------------------------------------
for j=1:nvar
    Yvec = Y(:,j);
    if ML == 0 %% Analytical ML
    OLSout = OLSmodel(Yvec,X,0);
    aux = ['eq' num2str(j)];
    eval( ['FAVAR.' aux '.beta  = OLSout.beta;'] );  % bhats
    eval( ['FAVAR.' aux '.tstat = OLSout.tstat;'] ); % t-stats
    % compute t-probs
    tstat = zeros(ncoeff,1);
    tstat = OLSout.tstat;
    tout = tdis_prb(tstat,nobse-ncoeff);
    eval( ['FAVAR.' aux '.tprob = tout;'] );         % t-probs
    eval( ['FAVAR.' aux '.resid = OLSout.resid;'] ); % resids 
    eval( ['FAVAR.' aux '.yhat  = OLSout.yhat;'] );  % yhats
    eval( ['FAVAR.' aux '.y     = Yvec;'] );         % actual y
    eval( ['FAVAR.' aux '.rsqr  = OLSout.rsqr;'] );  % r-squared
    eval( ['FAVAR.' aux '.rbar  = OLSout.rbar;'] );  % r-adjusted
    eval( ['FAVAR.' aux '.sige  = OLSout.sige;'] );  % standard error
    else %% Numerical ML
    [Meanmodel, volatilitymodel] = FAVAR_Univ_GARCH(Yvec,X);
    aux = ['eq' num2str(j)];
    eval( ['FAVAR.' aux '.beta  = Meanmodel.betas;'] );  % bhats
    Ft(:,j)  = Meanmodel.betas';
    Residuals(:,j)  = volatilitymodel.standardizedresidual;
    Volat(:,j)      = volatilitymodel.filteredvolatility;
    eval( ['FAVAR.' aux '.sige  = Meanmodel.sige;'] ); % standard errors
    eval( ['FAVAR.' aux '.yhat  = Meanmodel.yhat;'] );  % filtered mean
    eval( ['FAVAR.' aux '.residuals = Meanmodel.resid;'] ); % Residuals
    eval( ['FAVAR.' aux '.LLH   = Meanmodel.LLH;'] ); % Total LLH
    eval( ['FAVAR.' aux '.filteredvolatility   = volatilitymodel.filteredvolatility;']); % Total LLH
    eval( ['FAVAR.' aux '.standardizedresidual   = volatilitymodel.standardizedresidual;']);
    eval( ['FAVAR.' aux '.GARCHLLH   = volatilitymodel.LLH;']);
    eval( ['FAVAR.' aux '.GARCHpars  = volatilitymodel.pars;']);
    end 
end
%% Compute the matrix of coefficients & VCV
%--------------------------------------------------------------------------
FAVAR.GARCHresiduals = ML;
if ML == 0
    Ft = (X'*X)\(X'*Y);
    FAVAR.residuals = Y - X*Ft;
else
    FAVAR.residualsNonFilter = Y - X*Ft;
    FAVAR.residuals          = Residuals;
    FAVAR.volatility         = Volat;
end
FAVAR.Ft = Ft;
%SIGMA = (1/(nobse-ntotcoeff))*FAVAR.residuals'*FAVAR.residuals; % adjusted for # of estimated coeff per equation
SIGMA  = cov(FAVAR.residuals);
SSE  = sum((Y - X*Ft).^2);
FAVAR.sigma = SIGMA;
FAVAR.LLH = -(nobs/2)* (nvar*(1+log(2*pi)) + log(det(SIGMA)));
FAVAR.fit = X*Ft;
FAVAR.X = X;
FAVAR.Y = Y;
if nvar_ex > 0
    FAVAR.X_EX = X_EX;
end
if ML == 1
    if DCC == 1
    DCC = Estimate_DCC(FAVAR);
    FAVAR.DCC = DCC;
    else
        DCC = [];
    end
end
%% Companion matrix of Ft' and max eigenvalue
%--------------------------------------------------------------------------
F = Ft';
Fcomp = [F(:,1+const:nvar*nlag+const); eye(nvar*(nlag-1)) zeros(nvar*(nlag-1),nvar)];
FAVAR.Fcomp = Fcomp;
FAVAR.maxEig = max(abs(eig(Fcomp)));
FAVAR.Eigs = eig(Fcomp);

%% Initialize other results
%--------------------------------------------------------------------------
FAVAR.invA = [];  % inverse of the A matrix (need identification: see VARir/VARfevd)
FAVAR.S    = [];  % Orthonormal matrix (need identification: see SR)