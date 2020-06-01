clear all; clc
addpath(genpath('Data'));
addpath(genpath('Utils'));
clear all; clc;
set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
Information_Crit = 0;
MLGARCH          = 1;
HD               = 0;
Cholesky         = 1;
FFR1d            = 1;  % Switch to determine whether FFR should be assumed I(1) or levels I(0)
type             = 2;  % Switch 1 = PCA extraction of factors, 2 = Kalman Smoothing of factors. NB. Kalman is computationally very heavy.
det              = 0;
faclag           = 3;
nsteps           = 24; % Number of steps IRF
ConfLevel        = 68; % Confidence level percentage for bootstrap
if Cholesky == 1
    Order = [1,2,3];
else
    Order = [2,3,1];
end
var_numbers = [1, 19, 24, 66, 76, 80, 84,  104, 128, 129, 130 131];
              %1 = RPI, 19 = Capacity util, 24 = unrate, 66 = M2Real, 76 =
              %div yield, 80 = Tb3m, 84 = Tb10, 104 = Oil Price, 128 = Vola
var_names   = {'Real Personal Income', 'Capacity Utilization: Manufacturing',...
    'Civilian Unemployment Rate','Real M2 Money Stock','S&P Dividend Yield',...
    '3-Month Treasury Bill','10-Year Treasury Bond','Crude Oil Price','Volatility Index'}';
%% Full Sample 1959:2019
%-----------------------------Load data------------------------------------
PrepdataFull; 
if FFR1d == 1
    firstdstr = '1d';
else
    firstdstr = 'Level';
    ydata_st  = ydata_st(3:end,:);
end
Y = ydata_st;
X = x2; X_st = x2_st;
%%
%-----------------------Information Criterion------------------------------
if Information_Crit
    MaxFac = 10; MaxLag = 12; 
    [nfact, ~, ~]  = NumFacsDFM_FAVAR(X, MaxFac);
    saveas(gcf,['../Figures/FullNumFacs.pdf']);
    K = nfact(:,2500); %% 3000'th entry seems stable, K = number of factors;
    [NumLags, IC_table] = DetermineLagLengthFAVAR(K,faclag,MaxLag,Y,X_st,slowcode,det,type);
    plag = NumLags(3);         % plag is number of lags in the VAR part
else
    K = 3;
    plag = 3;
end
[FAVAR, FAVARopt] = EstimateReducedFAVAR(X, Y, faclag, plag, K, det, slowcode, type, namesXY,MLGARCH,Cholesky);
for i = 1:K
   vnames{i} = strcat("Factor ", char(string(i)));
end
vnames = [vnames';namesY];
if Cholesky == 0
FAVARopt.vnames = vnames([1,4,5,2,3,6]);
else
    FAVARopt.vnames = vnames;
end
PlotEigs(FAVAR,'EigsFullSample')
RMSEFull    = sqrt(mean(FAVAR.residuals.^2));
SummFullTab = GenerateSummaryStats(ydata,namesY);
if Cholesky == 1 % Ensure i dont by mistake overwrite my figure with a differnt ordering
FAVAR_Cusum(FAVAR,DATES,namesXY,K,['Full_Cusum',firstdstr],Cholesky);
end
PlotFilteredVolatilityFAVAR(FAVAR, FAVARopt, DATES, ['Filtered_Vola_full',firstdstr]);
Plot_FAVARfilter(FAVAR, FAVARopt, DATES, ['Filtered_FAVAR_Fullsample',firstdstr,'Lags',num2str(plag)]);
disp(sqrt(mean(FAVAR.residuals.^2)));
tests = ADFKPSS_Data(ydata_st, namesY, 1:12);
%% Structural Analysis
if Cholesky == 1
%---------------------------Cholesky SVAR----------------------------------
FAVARcholopt = GenerateOptionsFAVAR(1, ConfLevel, nsteps,vnames);
%-----------------------Impulse Response Functions-------------------------
[IRF_Chol, FAVAR_Chol] = VARir(FAVAR,FAVARcholopt);
% Compute error bands
[IRFINF_Chol,IRFSUP_Chol,IRFMED_Chol] = VARirband(FAVAR_Chol,FAVARcholopt);
save("IRFMedCholFull", "IRFMED_Chol")
% Plot VAR-responses
nlagstr = num2str(FAVAR_Chol.nlag);
name  = ['Chol',nlagstr,'Lags19592019',firstdstr];
%%
Modelname = ['3lagchol1959',firstdstr];
FAVARopt.pick = 1;
%set(0,'DefaultAxesColorOrder','remove');
%HD = VARhd(FAVAR_Chol);
%VARhdplot(HD,FAVARcholopt, FAVAR_Chol, DATES, Modelname);
else
%% Blanchard Quah Long Run Restrictions
%---------------------------Long Run SVAR----------------------------------
FAVARbqopt = GenerateOptionsFAVAR(2, ConfLevel, nsteps,vnames);
FAVARbqopt.standardized = 1;
FAVARbqopt.vnames = [{'Factor 1'}, {'IP'}, {'FFR'}, {'Factor 2'}, {'Factor 3'}, {'CPI'}]';
FAVARbqopt.snames = ["Factor 1", "Aggregate Supply", "Spending", "Factor 2", "Factor 3", "Expansionary monetary policy"]';
%-----------------------Impulse Response Functions-------------------------
[IRF_bq, FAVAR_bq] = VARir(FAVAR,FAVARbqopt);
% Compute error bands
[IRFINF_bq,IRFSUP_bq,IRFMED_bq] = FAVARirband(FAVAR_bq,FAVARbqopt);
save("IRFMedbqFull", "IRFMED_bq");
% Plot VAR-responses
nlagstr = num2str(FAVAR_bq.nlag);
name  = ['bq',nlagstr,'Lags19592019',firstdstr];
FAVAR_Resp(FAVAR_bq, FAVARbqopt, IRFMED_bq, IRFINF_bq, IRFSUP_bq,name, FFR1d,false);
%% Historical Decomposition
if HD == 1;
Modelname = ['3lagbq1959',firstdstr];
FAVARbqopt.pick = 1;
set(0,'DefaultAxesColorOrder','remove');
HD = VARhd(FAVAR_bq);
VARhdplot(HD,FAVARbqopt, FAVAR_bq, DATES, Modelname);
end
end