clearvars -except nrun; clc
addpath(genpath('Data'));
addpath(genpath('Utils'));
addpath(genpath('SVARToolbox'));
set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
if exist('nrun','var') == 0 % Checks if this is the first run;
   nrun = 0; 
end
%%
nrun             = 0;  % Set foor speed;
GARCH            = 1;  % switch GARCH(1,1) Residuals on/off; very computationally demanding
Cholesky         = 1;  % 1: Cholesky or 0: Blanchard/Quah
HD               = 0;  % Switch historical Decomposition
VARplot          = 0;  % Plot a comparison between FAVAR and VAR
FFR1d            = 0;  % Switch to determine whether FFR should be assumed I(1) or levels I(0)
type             = 2;  % Switch 1 = PCA extraction of factors, 2 = Kalman Smoothing of factors. NB. Kalman is computationally heavy.
plag             = 13;  % Number of lags in the VAR
faclag           = 3;  % Number of lags in the DFM (Measurement eq.)
K                = 3;  % Number of factors extracted
det              = 0;  % No constant
nsteps           = 48; % Number of steps IRF

if Cholesky == 1
    Order            = [1,2,3]; % 1: Inflation, 2: Industrial Prod, 3: FFR
                                % Cholesky order [1,2,3];
                                % BQ order: [2,3,1];
else
    Order            = [2,3,1];
end
ConfLevel        = 68; % Confidence level percentage for bootstrap
var_numbers = [1, 18, 23, 65, 75, 78, 82,  102, 125, 126, 127, 128];
              %1 = RPI, 18 = Capacity util, 23 = unrate, 65 = M2Real, 75 =
              %div yield, 78 = Tb3m, 82 = Tb10, 102 = Oil Price, 125 = Vola
var_names   = {'Real Personal Income', 'Capacity Utilization: Manufacturing',...
    'Civilian Unemployment Rate','Real M2 Money Stock','S&P Dividend Yield',...
    '3-Month Treasury Bill','10-Year Treasury Bond','Crude Oil Price','Volatility Index'}';
%% Cholesky ID Bernanke sample
%-----------------------------Load data------------------------------------
PrepDataBernanke;
%TableB = GenerateSummaryStats(ydata, namesXY(end-2:end,1));
%table2latex(TableB)
%%
var_names  = [var_names; namesY];
if FFR1d == 1
    if Order == [1,2,3]
    ydata(:,3) = [NaN; ydata(2:end,3)-ydata(1:end-1,3)];
    %ydata_st(:,3) = [NaN; ydata_st(2:end,3)-ydata_st(1:end-1,3)];
    %ydata = ydata(3:end);
    ydata_st = (ydata(3:end,:)-mean(ydata(3:end,:)))./std(ydata(3:end,:));
    else
    ydata(:,2) = [NaN; ydata(2:end,2)-ydata(1:end-1,2)];
    ydata_st = (ydata(3:end,:)-mean(ydata(3:end,:)))./std(ydata(3:end,:));
    %ydata_st(:,2) = (ydata_st(3:end,2)-mean(ydata_st(3:end,2)))./std(ydata_st(3:end,2));
    ydata = ydata(3:end,:);
    %ydata_st = ydata_st(3:end,:);
    end
    Y = ydata_st;
    firstdstr = '1d';
    DATES = DATES(3:end,:);
else
    Y = ydata_st(3:end,:);
    %Y = ydata(3:end,:);
    DATES = DATES(3:end,:);
    firstdstr = 'Level';
end
%%
X = x2_st;
clearvars vnames
for i = 1:K
   vnames{i} = strcat('Factor ', char(string(i)));
end
vnames = [vnames';namesY];
[FAVAR, FAVARopt] = EstimateReducedFAVAR(X, Y, faclag, plag, K, det, slowcode, type, namesXY, GARCH,Cholesky);
FAVAR_Cusum(FAVAR,DATES,namesXY,K,['Bernanke_Cusum',firstdstr],Cholesky);
FAVARopt.vnames = vnames;
%%
if Cholesky == 1;
[FAVAR, FAVARopt] = EstimateReducedFAVAR(X, Y, faclag, plag, K, det, slowcode, type, namesXY,GARCH,Cholesky);
FAVARopt.nsteps = 48;
FAVARopt.ndraws = 1000;
FAVARopt.pctg = 68;
FAVARopt.NegativeShockIndex = [];
%% Coefficient of determinination Factors, PCA vs. Kalman DFM
% Kalman
[FAVAR, FAVARopt] = EstimateReducedFAVAR(X, Y, faclag, plag, K, det, slowcode, type, namesXY,GARCH,Cholesky);
% PCA
[FAVAR_PCA, FAVAR_PCopt] = EstimateReducedFAVAR(X, Y, 1, plag, K, det, slowcode, 1, namesXY,0,Cholesky);
FAVARopt.vnames = vnames;
PCA_V_SmootherTab = table(FAVAR_PCA.DFM.R2s.Variables, FAVAR.DFM.R2s.Variables);
PCA_V_SmootherTab.Properties.VariableNames = ["PCA", "Kalman Smoothed"];
PCA_V_SmootherTab = PCA_V_SmootherTab(var_numbers,1:2);
Difference =  FAVAR.DFM.R2s.Variables - FAVAR_PCA.DFM.R2s.Variables;
mean(Difference)
R2s = (FAVAR.DFM.R2s.Variables-FAVAR_PCA.DFM.R2s.Variables);
R2s(R2s==0) = NaN; 
mean(R2s, 'omitnan');
%%
tic
FAVARCholopt = GenerateOptionsFAVAR(1, ConfLevel, nsteps,vnames);
FAVARCholopt.ndraws = 1000;
[FAVAR, FAVARopt] = EstimateReducedFAVAR(X, Y, faclag, plag, K, det, slowcode, type, namesXY,GARCH,Cholesky);
% Compute IRF
[IRF_Chol, FAVAR_Chol] = VARir(FAVAR,FAVARCholopt);
% Compute error bands
[IRFINF_Chol,IRFSUP_Chol,IRFMED_Chol, INFx, SUPx, MEDx] = FAVARirband(FAVAR_Chol,FAVARCholopt);
toc
save('IRFMedCholBern', 'IRFMED_Chol')
% Plot 
nlagstr = num2str(FAVAR_Chol.nlag);
name  = ['chol',nlagstr,'LagsBernankeSample',firstdstr];
FAVAR_Resp(FAVAR_Chol, FAVARCholopt, IRFMED_Chol, IRFINF_Chol, IRFSUP_Chol,name, FFR1d, false);
%%
if plag == 3 & K == 3
figure;
subplot(2,1,1)
yt = plot(exp(cumsum(MEDx(:,find(strcmp(namesXY, 'S&P 500')),6)))-1);
hold on
title([{'S\&P Common Stock Price Index: Composite';'to monetary policy'}]);ylabel({'Deviation from mean, std. units'});
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
              xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
xlim([1, nsteps]); yline(0,'--');
xt = fill([1:nsteps, fliplr(1:nsteps)], [exp(cumsum(INFx(:,find(strcmp(namesXY, 'S&P 500')),6)))'-1, fliplr(exp(cumsum(SUPx(:,find(strcmp(namesXY, 'S&P 500')),6))')-1)],'k','FaceAlpha',0.05,'HandleVisibility','off');
lg = legend([yt, xt], 'Median Response','Confidence interval $\pm 1 \sigma$','interpreter','latex','FontSize',14,'location','best');
subplot(2,1,2)
plot(cumsum(MEDx(:,find(strcmp(namesXY, 'S&P div yield')),6)));
hold on
fill([1:nsteps, fliplr(1:nsteps)], [cumsum(INFx(:,find(strcmp(namesXY, 'S&P div yield')),6))', fliplr(cumsum(SUPx(:,find(strcmp(namesXY, 'S&P div yield')),6))')],'k','FaceAlpha',0.05,'HandleVisibility','off');
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
              xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
title([{'S\&P Common Stock Dividend Yield';'to monetary policy'}]);ylabel({'Deviation from mean, std. units'});
xlim([1, nsteps]); yline(0,'--');xlabel('Periods','interpreter','latex','FontSize',14); 
xlabel('Periods','interpreter','latex','FontSize',14);
name = ['IRFsDFMSP5003LagsBern',firstdstr];
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
saveas(gcf, ['../Figures/',name,'.pdf'])
end
%%
prefix = ['CholB', firstdstr];
IRF_DFM(FAVAR, FAVARCholopt, MEDx, INFx, SUPx,  tcode, var_numbers, var_names, prefix, Order, FFR1d);
%% Historical Decomposition
if HD == 1;
Modelname = ['3lagchol19592001',firstdstr];
FAVARCholopt.pick = 1;
set(0,'DefaultAxesColorOrder','remove');
HD = VARhd(FAVAR_Chol);
VARhdplot(HD,FAVARCholopt, FAVAR_Chol, DATES, Modelname);
end
else
%% BQ
%---------------------------Long Run SVAR----------------------------------
Conf = 68;
nsteps = 48;
Y = ydata_st;
% Set options some options for IRF calculation and bands
[FAVAR, FAVARopt] = EstimateReducedFAVAR(X, Y, faclag, plag, K, det, slowcode, type, namesXY,GARCH,Cholesky, DCC);
FAVARbqopt = GenerateOptionsFAVAR(2, Conf, nsteps,vnames);
FAVARbqopt.standardized = 0;
% Compute IRF
[IRF_bq, FAVAR_bq] = VARir(FAVAR,FAVARbqopt);
% Compute error bands
[IRFINF_bq,IRFSUP_bq,IRFMED_bq,INFx,SUPx,MEDx] = FAVARirband(FAVAR,FAVARbqopt);
if K == 3 || plag == 3
    save("IRFMedbqBern", "IRFMED_bq")
end
% Plot 
nlagstr = num2str(FAVAR.nlag);
name  = ['bq',nlagstr,'LagsBernankeSample',firstdstr];
FAVAR_Resp(FAVAR_bq, FAVARbqopt, IRFMED_bq, IRFINF_bq, IRFSUP_bq,name, FFR1d, 0);
%% DFM IRF
prefix = ['BQb', firstdstr];
IRF_DFM(FAVAR, FAVARbqopt, MEDx, INFx, SUPx,  tcode, var_numbers, var_names, prefix, Order, FFR1d)
%% Historical Decomposition
if HD == 1;
Modelname = ['3lagbq19592001',firstdstr];
FAVARbqopt.pick = 1;
set(0,'DefaultAxesColorOrder','remove');
HD = VARhd(FAVAR_bq);
VARhdplot(HD,FAVARbqopt, FAVAR_bq, DATES, Modelname);
end
end
%% Different specifications:
% This section estimates the FAVAR with PCA and DFM and plots responses:
% p = 3, p = 13
if Cholesky == 1;
    if nrun == 0
Conf = 68;
[VAR3, VARopt] = VARmodel(Y,3,det,1); 
[VAR13, ~] = VARmodel(Y,13,det,1); 
VARopt.nsteps = 48;
VARopt.ndraws = 1000;         
[FAVAR3lag, ~] = EstimateReducedFAVAR(X, Y, faclag, 3, K, det, slowcode, 2, namesXY,1,1);
[FAVAR13lag, ~] = EstimateReducedFAVAR(X, Y, faclag, 13, K, det, slowcode, 2, namesXY,1,1);
[FAVARPCA3lag, ~] = EstimateReducedFAVAR(X, Y, 1, 3, K, det, slowcode, 1, namesXY,1,1);
[FAVARPCA13lag, ~] = EstimateReducedFAVAR(X, Y, 1, 13, K, det, slowcode, 1, namesXY,1,1);
FAVARDiffopt = GenerateOptionsFAVAR(1, Conf, nsteps,vnames);
FAVARDiffopt.NegativeShockIndex = [];
[IRFINF_3l,IRFSUP_3l,IRFMED_3l,~,~,~] = FAVARirband(FAVAR3lag,FAVARDiffopt);
[IRFINF_13l,IRFSUP_13l,IRFMED_13l,~,~,~] = FAVARirband(FAVAR13lag,FAVARDiffopt);
[IRFINF_PCA3l,IRFSUP_PCA3l,IRFMED_PCA3l,~,~,~] = FAVARirband(FAVARPCA3lag,FAVARDiffopt);
[IRFINF_PCA13l,IRFSUP_PCA13l,IRFMED_PCA13l,~,~,~] = FAVARirband(FAVARPCA13lag,FAVARDiffopt);
VARopt.NegativeShockIndex = [];
[~,~,IRFMED_VAR3] = VARirband(VAR3,VARopt);
[~,~,IRFMED_VAR13] = VARirband(VAR13,VARopt);
t = 1:48;
CPIResp = exp(cumsum([IRFMED_3l(:,4,6), IRFMED_13l(:,4,6), IRFMED_PCA3l(:,4,6), IRFMED_PCA13l(:,4,6), IRFMED_VAR3(:,1,3), IRFMED_VAR13(:,1,3)]))-1;
IPResp  = exp(cumsum([IRFMED_3l(:,5,6), IRFMED_13l(:,5,6), IRFMED_PCA3l(:,5,6), IRFMED_PCA13l(:,5,6), IRFMED_VAR3(:,2,3), IRFMED_VAR13(:,2,3)]))-1;
if FFR1d == 1
    FFRResp = cumsum([IRFMED_3l(:,6,6), IRFMED_13l(:,6,6), IRFMED_PCA3l(:,6,6), IRFMED_PCA13l(:,6,6), IRFMED_VAR3(:,3,3), IRFMED_VAR13(:,3,3)]);
else
FFRResp = [IRFMED_3l(:,6,6), IRFMED_13l(:,6,6), IRFMED_PCA3l(:,6,6), IRFMED_PCA13l(:,6,6), IRFMED_VAR3(:,3,3), IRFMED_VAR13(:,3,3)];
end
%%
figure
set(0,'DefaultAxesColorOrder',[0 0 0]);
subplot(1,3,1);
plot(t,CPIResp(t,1),'-',t,CPIResp(t,2),'--',t,CPIResp(t,3),':',t,CPIResp(t,4),'-.','LineWidth',1.5,'MarkerSize',3);
hold on
plot(t,CPIResp(t,5),'-o',t,CPIResp(t,6),'-*','Color',uint8([50,50,50]));
yline(0,'-','LineWidth',1);
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
              xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
%legend(["DFM, $p=3$","DFM, $p = 13$","PCA, $p=3$","PCA, $p=13$"],'interpreter','latex','location','best')
title('Monetary Policy shock to CPI');
subplot(1,3,2);
plot(t,IPResp(t,1),'-',t,IPResp(t,2),'--',t,IPResp(t,3),':',t,IPResp(t,4),'-.','LineWidth',1.5,'MarkerSize',3);
hold on
plot(t,IPResp(t,5),'-o',t,IPResp(t,6),'-*','Color',uint8([50,50,50]));
yline(0,'-','LineWidth',1);
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
              xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
%legend(["DFM, $p=3$","DFM, $p = 13$","PCA, $p=3$","PCA, $p=13$"],'interpreter','latex','location','best')
title('Monetary Policy shock to IP')
subplot(1,3,3);
plot(t,FFRResp(t,1),'-',t,FFRResp(t,2),'--',t,FFRResp(t,3),':',t,FFRResp(t,4),'-.','LineWidth',1.5,'MarkerSize',3);
hold on
plot(t,FFRResp(t,5),'-o',t,FFRResp(t,6),'-*','Color',uint8([50,50,50]));
title('Monetary Policy shock to FFR');
yline(0,'-','LineWidth',1);
              yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
              xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
set(gcf, 'unit', 'inches');
lg = legend(["DFM, $p=3$","DFM, $p = 13$","PCA, $p=3$","PCA, $p=13$", "VAR, $p=3$", "VAR, $p=13$"],'interpreter','latex','Orientation','horizontal');
set(lg,'Units','inches',...
    'Position',[3.5 0.1 1.0 0.23],...
    'Interpreter','latex',...
    'FontSize',14);
set(gcf, 'PaperPosition', [0 0 40 18]);
set(gcf, 'PaperSize', [40 18]);
saveas(gcf,['../Figures/MonetaryPolicyDiffModels',firstdstr,'.pdf']);
end
nrun = 1;
end
