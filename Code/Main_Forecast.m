% Script conducting forecasts in the FAVAR and compares to benchmark
% forecasts. A series of CW-tests are produces as well as relative RMSE are
% calculated and saved in a table. Additional out-of-sample measures are
% calculated including sign-accuracy and estimated FC bias.
% All forecasts are generated in a rolling window style.

% Rasmus M. Jensen 2020
clear all;
tic;
addpath(genpath('Data'));
addpath(genpath('Utils')); clear all; clc;
set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
%% Settings
tic
Order = [1,2,3];          % Order of variables 1: infl, 2: IP, 3 FFR.
FFR1d = 1;                % FFR: I(0) I(1) 
horizons = [1, 3, 6, 12]; % Horizons to consider
NfacPCA = 3; % Number of factors PCA.
plag    = 3; % Number of lags in the VAR (state-equation)
faclag  = 3; % Number of factor lags in the DFM (signal equation)
K       = 3; % Number of factors DFM
det     = 0; % Deterministics, 0 none, 1 constant
%%
PrepdataFull;
YFull = ydata_st;
XFull = x2_st;
DATEY = DATES;
WinSize = 120;  
NamesXYfull= namesXY;
type    = 2; % dont change.
%% FullSample
clearvars FEn Fitn Obsn Daten
iter = 0;
    for i=WinSize+2:size(YFull,1)-(horizons(end)+1) % Loop over sample
        iter = iter + 1;
        [FEn(iter,:,:), Fitn(iter,:,:), Obsn(iter,:,:),...
            Daten(iter,:)] = RollingWindowFAVARComp(XFull, YFull,faclag, plag,...
            K, det, slowcode,DATEY,WinSize,i,type,horizons,false,NamesXYfull);
        [FE1Fn(iter,:,:), Fit1Fn(iter,:,:), Obs1Fn(iter,:,:),...
            ~] = RollingWindowFAVARComp(XFull, YFull,faclag, plag,...
            1, det, slowcode,DATEY,WinSize,i,type,horizons,false,NamesXYfull);
        [FE1l(iter,:,:), Fit1l(iter,:,:), Obs1l(iter,:,:),...
            ~] = RollingWindowFAVARComp(XFull, YFull,faclag, 1,...
            K, det, slowcode,DATEY,WinSize,i,type,horizons,false,NamesXYfull);
    end
clearvars FEPCAn FitPCAn ObsPCAn DatePCAn
iter = 0;
    for i=WinSize+2:size(YFull,1)-(horizons(end)+1) % Loop over sample
        iter = iter + 1;
        [FEPCAn(iter,:,:), FitPCAn(iter,:,:), ObsPCAn(iter,:,:),...
            DatePCAn(iter,:)] = RollingWindowFAVARComp(XFull, YFull,1, plag,...
            NfacPCA, det, slowcode,DATEY,WinSize,i,1,horizons,false,NamesXYfull);
        [FE8PCAn(iter,:,:), Fit8PCAn(iter,:,:), Obs8PCAn(iter,:,:),...
            Date8PCAn(iter,:)] = RollingWindowFAVARComp(XFull, YFull,1, plag,...
            8, det, slowcode,DATEY,WinSize,i,1,horizons,false,NamesXYfull);
    end
%%
% AR(1), RW, MM, VAR
clearvars Ytp1MM phi FEAR1 Ytp1AR FEMM Ytp1RW FERW VARFE VARPred VARObs Vardates
Y = YFull;
DATES = DATEY;
for k = 1:3
    hiter = 0;
    for h = horizons
        hiter = hiter + 1;
        iter = 0;
        for t= WinSize + 2:size(Y,1)-(horizons(end)+1)
        iter = iter + 1;
        %[Ytp1AR(k,iter,hiter), phi(k,iter,hiter), FEAR1(k,iter,hiter), OBSAR1(k,iter,hiter)] = AR_FC(Y(:,k), t,WinSize,h);
        [Ytp1AR(k,iter,hiter), ~, FEAR1(k,iter,hiter), OBSAR1(k,iter,hiter)] = AR_FC(Y(:,k), t,WinSize,h);
        [Ytp1MM(k,iter,hiter), FEMM(k,iter,hiter)] = MeanModelFC(Y(:,k), t, WinSize, h);
        [Ytp1RW(k,iter,hiter), FERW(k,iter,hiter), OBSrw(k,iter,hiter)] = RandomwalkForecast(Y(:,k), t, h);
        end
    end
end
iter = 0;
for t= WinSize + 2:size(Y,1)-(horizons(end)+1)
    iter = iter + 1;
    [VARFE(:,:,iter), VARPred(:,:,iter), VARObs(:,:,iter), Vardates(iter,:)] = FC_VAR(YFull, plag, t, WinSize, DATES,horizons);
end
%% Clark West test
Modelnames = {'FAVAR, PCA k=3', 'FAVAR 1 lag', 'FAVAR 1 Factor','AR1','VAR','Random Walk'};
for h = 1:4
    for series = 1:3 
      Predall   = [squeeze(FitPCAn(:,series,h)),squeeze(Fit1l(:,series,h)),squeeze(Fit1Fn(:,series,h)),squeeze(Ytp1AR(series,:,h))',squeeze(VARPred(series,h,:)),squeeze(Ytp1RW(series,:,h))'];
      PredBench = squeeze(Fitn(:,series,h)); % squeeze(Ytp1RW(series,:,h))';%;squeeze(Ytp1AR(series,:,h))';
      OBS       = squeeze(Obsn(:,series,h));
        for model = 1:6
            [cw_stat(model,series,h),cw_pval(model,series,h)] = cw(Predall(:,model), PredBench,OBS);
        end
    end
end
cwtab = splitvars(table(round(reshape(cw_pval,6,12),3)));
cwtab.Properties.RowNames = Modelnames'; 
disp(cwtab); % Disp. pvalues of the Clark & West (2007) test
%table2latex(cwtab, '../Tables/CW_tab.tex');
%%
FEn(FEn==0) = NaN;
FEPCAn(FEPCAn==0) = NaN;
FE8PCAn(FE8PCAn==0) = NaN;
FEAR1(FEAR1==0) = NaN; 
FEMM(FEMM==0) = NaN; 
FERW(FERW==0) = NaN; 
VARFE(VARFE==0) = NaN; 
FE1Fn(FE1Fn==0) = NaN;

RMSFEtabPCA8     = splitvars(table(squeeze(sqrt(mean(FE8PCAn.^2,1,'omitNaN'))))); 
RMSFEtabAR1     = splitvars(table(squeeze(sqrt(mean(FEAR1.^2,2,'omitNaN'))))); 
RMSFEtabVAR     = splitvars(table(squeeze(sqrt(mean(VARFE.^2,3,'omitNaN'))))); 
RMSFEtabMM      = splitvars(table(squeeze(sqrt(mean(FEMM.^2,2,'omitNaN')))));  
RMSFEtabRW      = splitvars(table(squeeze(sqrt(mean(FERW.^2,2,'omitNaN'))))); 
RMSFEtabPCA1959 = splitvars(table(squeeze(sqrt(mean(FEPCAn.^2,'omitNaN'))))); % Mean squared forecast error of the FAVAR
RMSFEtab1959    = splitvars(table(squeeze(sqrt(mean(FEn.^2,'omitNaN'))))); % Mean squared forecast error of the FAVAR
RMSEtabFAVAR1F  = splitvars(table(squeeze(sqrt(mean(FE1Fn.^2,'omitNaN'))))); % Mean squared forecast error of the FAVAR
RMSEtabFAVAR1L  = splitvars(table(squeeze(sqrt(mean(FE1l.^2,'omitNaN'))))); % Mean squared forecast error of the FAVAR
%% Full Measures
ARPerf        = RMSFEtabAR1.Variables./RMSFEtab1959.Variables;     % AR/RW
RWPerf        = RMSFEtabRW.Variables./RMSFEtab1959.Variables;      % RW/RW
VARPerf       = RMSFEtabVAR.Variables./RMSFEtab1959.Variables;     % VAR/RW
FAVARPerf     = RMSFEtab1959.Variables./RMSFEtab1959.Variables;    % FAVAR/RW
PCAFAVARPerf  = RMSFEtabPCA1959.Variables./RMSFEtab1959.Variables; % PCA-FAVAR/RW
PCA8FAVARPerf = RMSFEtabPCA8.Variables./RMSFEtab1959.Variables;    % FAVAR/RW
FAVARPerf1F   = RMSEtabFAVAR1F.Variables./RMSFEtab1959.Variables;  % FA(1)VAR/RW
FAVARPerf1L   = RMSEtabFAVAR1L.Variables./RMSFEtab1959.Variables;  % FA(1)VAR/RW

%%
InflFC = [RWPerf(1,:); FAVARPerf(1,:);FAVARPerf1L(1,:); FAVARPerf1F(1,:);PCA8FAVARPerf(1,:);PCAFAVARPerf(1,:);  VARPerf(1,:); ARPerf(1,:)];
GDPFC  = [RWPerf(2,:); FAVARPerf(2,:);FAVARPerf1L(2,:); FAVARPerf1F(2,:) ;PCA8FAVARPerf(2,:);PCAFAVARPerf(2,:); VARPerf(2,:); ARPerf(2,:)];
IntRFC = [RWPerf(3,:); FAVARPerf(3,:);FAVARPerf1F(3,:); FAVARPerf1F(3,:) ;PCA8FAVARPerf(3,:);PCAFAVARPerf(3,:); VARPerf(3,:); ARPerf(3,:)];
%% 
PerfTable = table(InflFC, GDPFC, IntRFC); PerfTable.Variables = round(PerfTable.Variables, 2);
PerfTable.Properties.RowNames = ["RW", "FAVAR-DFM 3 factors", "FAVAR-1-lag","FAVAR-DFM, 1 factor","FAVAR-8PCA", "FAVAR-PCA", "VAR", "AR"]
%table2latex(PerfTable,'../Tables/ForeCastMeas.tex')
%%
squeeze(mean(FEn,'omitnan'));
squeeze(mean(FERW,2, 'omitnan'));
squeeze(mean(FEAR1,2,'omitnan'));
squeeze(mean(FEPCAn,'omitnan'));
squeeze(mean(VARFE,3,'omitnan'));
%%
clearvars FAVAR1Fsign FAVARPCAsign FAVARsign RWsign VARsign
%%
for j = 1:size(OBSAR1, 1)
    for h = 1:size(OBSAR1, 3)
        for i = 1:size(OBSAR1, 2)
            if sign(OBSAR1(j,i,h)) == sign(Ytp1AR(j,i,h));
                ARsign(j,i,h) = 1;
            else
                ARsign(j,i,h) = 0;
            end
            if sign(Ytp1MM(j,i,h)) == sign(OBSAR1(j,i,h));
                MMsign(j,i,h) = 1;
            else
                MMsign(j,i,h) = 0;
            end
            if sign(Ytp1RW(j,i,h)) == sign(OBSrw(j,i,h));
                RWsign(j,i,h) = 1;
            else
                RWsign(j,i,h) = 0;
            end
        end
    end
end
for i = 1:size(Fitn,1)
    for j = 1:size(Fitn,2)
        for h = 1:size(Fitn,3)
            if sign(Fitn(i,j,h)) == sign(Obsn(i,j,h))
                FAVARsign(i,j,h) = 1;
            else
                FAVARsign(i,j,h) = 0;
            end
            if sign(Fit1Fn(i,j,h)) == sign(Obs1Fn(i,j,h))
                FAVAR1Fsign(i,j,h) = 1;
            else
                FAVAR1Fsign(i,j,h) = 0;
            end
            if sign(FitPCAn(i,j,h)) == sign(ObsPCAn(i,j,h))
                FAVARPCAsign(i,j,h) = 1;
            else
                FAVARPCAsign(i,j,h) = 0;
            end
            if sign(VARPred(j,h,i)) == sign(ObsPCAn(i,j,h))
                VARsign(i,j,h) = 1;
            else
                VARsign(i,j,h) = 0;
            end
            if sign(Fit8PCAn(i,j,h))  == sign(Obs8PCAn(i,j,h))
                PCA8(i,j,h) = 1;
            else
                PCA8(i,j,h) = 0;
            end
            if sign(Fit1l(i,j,h))  == sign(Obs1l(i,j,h))
                DFM1l(i,j,h) = 1;
            else
                DFM1l(i,j,h) = 0;
            end
        end
    end
end
%%
FAVARsign    = squeeze(mean(FAVARsign,'omitnan'));
FAVAR1Fsign  = squeeze(mean(FAVAR1Fsign,'omitnan'));
FAVARPCAsign = squeeze(mean(FAVARPCAsign,'omitnan'));
FAVAR1lsign  = squeeze(mean(DFM1l,'omitnan'));
FAVARPCA8sign = squeeze(mean(PCA8,'omitnan'));
MMsign       = squeeze(mean(MMsign,2,'omitnan'));
RWsign       = squeeze(mean(RWsign,2,'omitnan'));
ARSign       = squeeze(mean(ARsign,2,'omitnan'));
VARsign      = squeeze(mean(VARsign,1,'omitnan'));
    %%
INFsign = [FAVARsign(1,:); FAVAR1Fsign(1,:);FAVAR1lsign(1,:);FAVARPCA8sign(1,:); FAVARPCAsign(1,:);VARsign(1,:);ARSign(1,:);RWsign(1,:)];
Prodsign = [FAVARsign(2,:); FAVAR1Fsign(2,:);FAVAR1lsign(2,:);FAVARPCA8sign(2,:); FAVARPCAsign(2,:);VARsign(2,:);ARSign(2,:);RWsign(2,:)];
INTsign = [FAVARsign(3,:); FAVAR1Fsign(3,:);FAVAR1lsign(3,:);FAVARPCA8sign(3,:); FAVARPCAsign(3,:);VARsign(3,:);ARSign(3,:);RWsign(3,:)];
signtab = round([INFsign, Prodsign, INTsign], 3)*100;
signtab = splitvars(table(signtab));
signtab.Properties.RowNames = ["DFM, k=3", "DFM k=1","DFM  l=1","PCA k = 8", "PCA k = 3", "VAR", "AR", "RW"]
%table2latex(signtab,'../Tables/SignTable.tex')
toc