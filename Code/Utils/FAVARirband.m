function [INF,SUP,MED,INFx,SUPx,MEDx] = FAVARirband(VAR,VARopt)
% =======================================================================
% Calculate confidence intervals for impulse response functions computed
% with VARir - Modified to accomodate a DFM in the FAVAR framework
% --------------------------------Input------------------------------------
%   VAR: structure, result of VARmodel -> VARir function
%   VARopt: options of the VAR (see VARopt from VARmodel)
%---------------------------------Output-----------------------------------
%   INF(t,j,k): lower confidence band (t steps, j variable, k shock)
%   SUP(t,j,k): upper confidence band (t steps, j variable, k shock)
%   MED(t,j,k): median response (t steps, j variable, k shock)
% DFM-Responses:
%   INFx(t,j,k): lower confidence band (t steps, j variable, k shock)
%   SUPx(t,j,k): upper confidence band (t steps, j variable, k shock)
%   MEDx(t,j,k): median response (t steps, j variable, k shock)

% Ambrogio Cesa Bianchi, August 2015
% Modified: Rasmus M. Jensen 2020
%% Parse inputs
if ~exist('VAR','var')
    error('You need to provide VAR structure, result of VARmodel');
end
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
%% Retrieve and initialize variables 
nsteps            = VARopt.nsteps;
ndraws            = VARopt.ndraws;
pctg              = VARopt.pctg;
method            = VARopt.method;
Ft                = VAR.Ft;  % rows are coefficients, columns are equations
nvar              = VAR.nvar;
nvar_ex           = VAR.nvar_ex;
nlag              = VAR.nlag;
nlag_ex           = VAR.nlag_ex;
const             = VAR.const;
nobs              = VAR.nobs;
resid             = VAR.residuals;
ENDO              = VAR.ENDO;
EXOG              = VAR.EXOG;
L                 = VAR.DFM.L;
nStd              = VARopt.nStd;
CIMethod          = VARopt.CIMethod;
GARCHresiduals    = VAR.GARCHresiduals;

INF = zeros(nsteps,nvar,nvar);
SUP = zeros(nsteps,nvar,nvar);
MED = zeros(nsteps,nvar,nvar);
%BAR = zeros(nsteps,nvar,nvar); % Mean response removed

%% Create the matrices for the loop
%===============================================
y_artificial = zeros(nobs+nlag,nvar);


%% Loop over the number of draws
%===============================================

tt = 1; % numbers of accepted draws
ww = 1; % index for printing on screen
while tt<=ndraws
    
    % Display number of loops
    if tt==10*ww
        disp(['Loop ' num2str(tt) ' / ' num2str(ndraws) ' draws'])
        ww=ww+1;
    end

%% STEP 1: choose the method and generate the residuals
    if GARCHresiduals == 1
        u   = resid(ceil(size(resid,1)*rand(nobs,1)),:);
        sig = VAR.volatility;
        u   = u .* sig; % Reshuffles primitive residuals and corrects for vola.
    elseif strcmp(method,'bs')
    %if strcmp(method,'bs') 
        % Use the residuals to bootstrap: generate a random number bounded 
        % between 0 and # of residuals, then use the ceil function to select 
        % that row of the residuals (this is equivalent to sampling with replacement)
        u = resid(ceil(size(resid,1)*rand(nobs,1)),:);
    elseif strcmp(method,'wild')
        % Wild bootstrap based on simple distribution (~Rademacher)
        rr = 1-2*(rand(nobs,1)>0.5);
        u = resid.*(rr*ones(1,nvar));
    else
        error(['The method ' method ' is not available'])
    end

%% STEP 2: generate the artificial data
    % STEP 2.1: initial values for the artificial data
    % Intialize the first nlag observations with real data
    LAG=[];
    for jj = 1:nlag
        y_artificial(jj,:) = ENDO(jj,:);
        LAG = [y_artificial(jj,:) LAG]; 
    end
    % Initialize the artificial series and the LAGplus vector
    T = [1:nobs]';
    if const==0
        LAGplus = LAG;
    elseif const==1
        LAGplus = [1 LAG];
    elseif const==2
        LAGplus = [1 T(1) LAG]; 
    elseif const==3
        T = [1:nobs]';
        LAGplus = [1 T(1) T(1).^2 LAG];
    end
    if nvar_ex~=0
        LAGplus = [LAGplus VAR.X_EX(jj-nlag+1,:)];
    end
    % STEP 2.2: generate artificial series
    % From observation nlag+1 to nobs, compute the artificial data
    for jj = nlag+1:nobs+nlag
        for mm = 1:nvar
            % Compute the value for time=jj
            y_artificial(jj,mm) = LAGplus * Ft(1:end,mm) + u(jj-nlag,mm);
        end
        % now update the LAG matrix
        if jj<nobs+nlag
            LAG = [y_artificial(jj,:) LAG(1,1:(nlag-1)*nvar)];
            if const==0
                LAGplus = LAG;
            elseif const==1
                LAGplus = [1 LAG];
            elseif const==2
                LAGplus = [1 T(jj-nlag+1) LAG];
            elseif const==3
                LAGplus = [1 T(jj-nlag+1) T(jj-nlag+1).^2 LAG];
            end
            if nvar_ex~=0
                LAGplus = [LAGplus VAR.X_EX(jj-nlag+1,:)];
            end
        end
    end

%% STEP 3: estimate VAR on artificial data. 
    if nvar_ex~=0
        [VAR_draw, ~] = VARmodel(y_artificial,nlag,const,EXOG,nlag_ex);
    else
        [VAR_draw, ~] = VARmodel(y_artificial,nlag,const,0,0);
    end
%% STEP 4: calculate "ndraws" impulse responses and store them
    
    [irf_draw, VAR_draw] = VARir(VAR_draw,VARopt);  % uses options from VARopt, but companion etc. from VAR_draw
%% Step  4.5 Calculate the responses of either the dynamic factor model or the static factor model. Depending on whether PCA or DFM is used
        L = VAR.DFM.DynamicL; 
%%
if VAR_draw.maxEig<.9999
for i=1:VAR_draw.nvar
    if size(L,2) > size(VAR.DFM.L,2)
    foo = [irf_draw(1,:,i), zeros(1,size(L,2)-size(irf_draw,2))];
    else
    foo = irf_draw(1,:,i);
    end
    for ii = 1:size(irf_draw, 1)-1
            %IRFx(:,:,i,tt) = irf_draw(:,:,i) * L';
            %IRFxdraws(:,:,i,tt) = irf_draw(:,:,i) * L';
            
            IRFx(ii,:,i,tt)      = foo * L';
            IRFxdraws(:,:,i,tt)  = foo * L' ;
            if size(L,2) > size(VAR.DFM.L,2)
                foo = [irf_draw(1+ii,:,i), foo(1:nvar*2)];
            else
                foo = irf_draw(1+ii,:,i);
            end
            if ii == size(irf_draw, 1)-1;
                IRFx(ii+1,:,i,tt)      = foo * L';
            end
    end
end

        IRF(:,:,:,tt) = irf_draw;
        IRFdraws(:,:,:,tt)=irf_draw;
        tt=tt+1;
end
end
disp('-- Done!');
disp(' ');

%% Compute the error bands
if CIMethod == 1
MED(:,:,:) = median(IRF(:,:,:,:),4);
MEDx(:,:,:) = median(IRFx,4);
INF(:,:,:) = MED(:,:,:) + std(IRF,0,4) * (-nStd);
SUP(:,:,:) = MED(:,:,:) + std(IRF,0,4) * nStd;
INFx(:,:,:) = MEDx(:,:,:) + std(IRFx,0,4) * (-nStd);
SUPx(:,:,:) = MEDx(:,:,:) + std(IRFx,0,4) * nStd;
else
pctg_inf = (100-pctg)/2; 
pctg_sup = 100 - (100-pctg)/2;
INF(:,:,:)  = prctile(IRF(:,:,:,:),pctg_inf,4);
SUP(:,:,:)  = prctile(IRF(:,:,:,:),pctg_sup,4);
MED(:,:,:)  = prctile(IRF(:,:,:,:),50,4);
INFx(:,:,:) = prctile(IRFx(:,:,:,:),pctg_inf,4);
SUPx(:,:,:) = prctile(IRFx(:,:,:,:),pctg_sup,4);
MEDx(:,:,:) = prctile(IRFx(:,:,:,:),50,4);
end
%IRFs  = sort(IRF, 4); 
%IRFxs = sort(IRFx, 4);
%sizeD = size(IRFxs,4);
%INF(:,:,:) = IRFs(:,:,:,sizeD*.05);
%%INFx(:,:,:) = IRFxs(:,:,:,sizeD*.05);
%SUP(:,:,:) = IRFs(:,:,:,sizeD*.95);
%SUPx(:,:,:) = IRFxs(:,:,:,sizeD*.95);
%MED(:,:,:) = median(IRFs,4);
%MEDx(:,:,:) = median(IRFxs,4);
%BAR(:,:,:) = mean(IRF(:,:,:,:),4); Mean responses removed 2020, RMJ
%BARx(:,:,:) = mean(IRFx(:,:,:,:),4);
IRFdraws=IRF;
save IRFdraws IRFdraws
end


