function [FEVD, VAR, FAVFEVD] = FAVARfevd(VAR,VARopt)
% =======================================================================
% Compute FEVDs for a VAR model estimated with VARmodel. Three
% identification schemes can be specified: zero short-run restrictions,
% zero long run restrictions, and sign restrictions
% =======================================================================
% [FEVD, VAR] = VARfevd(VAR,VARopt)
%--------------------------------Input-------------------------------------
%   - VAR: structure, result of VARmodel function
%   - VARopt: options of the VAR (see VARopt from VARmodel)
%--------------------------------Output------------------------------------
%	- FEVD(t,j,k): matrix with 't' steps, the FEVD due to 'j' shock for 
%       'k' variable (Disabled)
%   - VAR: structure including VAR estimation results
%       * VAR.invA: identified contemporaneous A matrix
%   - FAVFEVD(t,j,n): array of t steps, the factor adjusted-FEVD for the
%   j'th shock to the in [x,y]_n variable.
% =======================================================================
% Based on: Ambrogio Cesa Bianchi, March 2015
%           ambrogio.cesabianchi@gmail.com

% Modified to accomodate the FAVAR-framework in 2020
% By Rasmus M. Jensen


%% Check inputs
%===============================================
if ~exist('VAR','var')
    error('You need to provide VAR structure, result of VARmodel');
end
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end


%% Retrieve and initialize variables 
%===============================================
nsteps = VARopt.nsteps;
ident  = VARopt.ident;

L     = VAR.DFM.L;
S     = VAR.S;
Fcomp = VAR.Fcomp;
nvar  = VAR.nvar;
sigma = VAR.sigma;

MSE   = zeros(nvar,nvar,nsteps);
MSE_j = zeros(nvar,nvar,nsteps);
PSI   = zeros(nvar,nvar,nsteps);
FEVD  = zeros(nsteps,nvar,nvar);
SE    = zeros(nsteps,nvar);


%% Compute the multipliers
%===============================================
VARjunk = VAR;
VARjunk.sigma = eye(nvar);
IRFjunk = VARir(VARjunk,VARopt);

% this loop is to save the multipliers for each period
for mm = 1:nvar
    PSI(:,mm,:) = reshape(IRFjunk(:,:,mm)',1,nvar,nsteps);
end


%% Calculate the contribution to the MSE for each shock (i.e, FEVD)
%===============================================

for mm = 1:nvar % loop for the shocks
    
    % Calculate Total Mean Squared Error
    MSE(:,:,1) = sigma;
    FAVMSE(:,:,1) = L * sigma * L';

    for kk = 2:nsteps
       MSE(:,:,kk) = MSE(:,:,kk-1) + PSI(:,:,kk)*sigma*PSI(:,:,kk)';
       FAVMSE(:,:,kk) = L * squeeze(MSE(:,:,kk)) * L';
    end
    
    % Get the matrix invA containing the structural impulses
    if strcmp(ident,'oir')
        [out, chol_flag] = chol(sigma);
        if chol_flag~=0; error('VCV is not positive definite'); end
        invA = out';
    elseif strcmp(ident,'bq')
        Finf_big = inv(eye(length(Fcomp))-Fcomp);   % from the companion
        Finf     = Finf_big(1:nvar,1:nvar);
        D        = chol(Finf*sigma*Finf')'; % identification: u2 has no effect on y1 in the long run
        invA = Finf\D;
    elseif strcmp(ident,'sr')
        [out, chol_flag] = chol(sigma);
        if chol_flag~=0; error('VCV is not positive definite'); end
        if isempty(S); error('Rotation matrix is not provided'); end
        invA = (out')*(S');
    end
    
    % Get the column of invA corresponding to the mm_th shock
    column = invA(:,mm);
    
    % Compute the mean square error
    MSE_j(:,:,1) = column*column';
    FAVMSE_J(:,:,1) = L * squeeze(MSE_j(:,:,1)) * L';
    for kk = 2:nsteps
        MSE_j(:,:,kk) = MSE_j(:,:,kk-1) + PSI(:,:,kk)*(column*column')*PSI(:,:,kk)';   
        FAVMSE_J(:,:,kk) = L * squeeze(MSE_j(:,:,kk)) * L';
    end

    % Compute the Forecast Error Covariance Decomposition
    FECD = MSE_j./MSE;
    FAVFECD = FAVMSE_J./FAVMSE;    
    
    % Select only the variance terms
    for nn = 1:nsteps
        for ii = 1:nvar
            FEVD(nn,mm,ii) = FECD(ii,ii,nn);
            SE(nn,:) = sqrt( diag(MSE(:,:,nn))' );
        end
        for kk = 1:size(L,1)
            FAVFEVD(nn,mm,kk) = FAVFECD(kk,kk,nn);
        end
    end
end

% Update VARopt
VAR.invA = invA;



