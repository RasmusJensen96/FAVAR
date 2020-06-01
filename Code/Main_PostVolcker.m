% FAVAR identification, assesment and structural analysis
% Script for the main sample. That is the post-Volcker period.
% Rasmus M. Jensen 2020
%%    
clear all;
tic;
addpath(genpath('Data'));
addpath(genpath('Utils'));
clear all; clc;
set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
set(0,'DefaultAxesColorOrder',[0 0 0]);

FREDhelp = readtable('totalDS.csv'); % Extended info table on the
           % factor panel, textmined from the McCracken & Ng (2015)-paper
           % in R
%% Settings
MLGARCH          = 1;         % Add GARCH(1,1) residuals to the model (DCC) (Maximum Liklihood)
Cholesky         = 1;         % 1 = Cholesky, else BQ.
HD               = 0;         % Historical Decompositions
FFR1d            = 1;         % Switch to determine whether FFR should be assumed I(1) or levels I(0)
Pre2009          = 0;         % Aux switch, 0 = 1959:2019, 1 = 1959:2008
type             = 2;         % Switch 1 = PCA extraction of factors, 2 = Kalman Smoothing of factors. NB. Kalman is computationally very heavy.
nsteps           = 24;        % Number of steps IRF
ConfLevel        = 68;        % Confidence level percentage for bootstrap
Information_Crit = 0;         % Switch 1: Lag length and numfacs by IC
faclag           = 3;         % Lags in the measurement equation
var_numbers = [1, 18, 23, 65, 75, 78, 82,  102, 125, 126, 127, 128];
              %1 = RPI, 18 = Capacity util, 23 = unrate, 65 = M2Real, 75 =
              %div yield, 78 = Tb3m, 82 = Tb10, 102 = Oil Price, 125 = Vola
var_names   = {'Real Personal Income', 'Capacity Utilization: Manufacturing',...
    'Civilian Unemployment Rate','Real M2 Money Stock','S\&P Dividend Yield',...
    '3-Month Treasury Bill','10-Year Treasury Bond','Crude Oil Price','Volatility Index'}';
if Cholesky == 1
    Order            = [1,2,3]; % 1: Inflation, 2: Industrial Prod, 3: FFR
                                % Cholesky order [1,2,3];
                                % BQ order: [2,3,1];
else
    Order            = [2,3,1];
end

%% Load 1983:2019 Data
%-----------------------------Load Data------------------------------------
% PrepData is a seperate script importing and preparing all dataseries,
% including removing outliers, standardizing and running the EM-algorithm
% it outputs all data neccesary for estimating the FAVAR, as well as
% an auxiliary vector of sample dates.
PrepData
TableA = GenerateSummaryStats(ydata, namesXY(end-2:end,1));
table2latex(TableA, './ResultsMATS/Summary1983')
%%
if FFR1d == 1
    firstdstr = '1d';
else
    firstdstr = 'Level';
end

if Pre2009 == 1
ydata_st = (ydata(1:300,:)-mean(ydata(1:300,:)))./std(ydata(1:300,:));
x2 = x2(1:300,:);
x2_st = (x2-mean(x2))./std(x2); 
DATES = DATES(1:300,:); DATESnm = DATESnm(1:300,:);
end

%------------------------Preliminary Stationarity check--------------------
lags = 1:12;
tests = ADFKPSS_Data(ydata_st, namesY, lags);
% ADF: 0 = unitroot, KPSS: 0 = trend-stationary
%% IC
%-----------------------Information Criterion------------------------------
Y = ydata_st; X_st = x2_st; X = x2; det = 0;
if Information_Crit
    MaxFac = 10; MaxLag = 12;
    [nfact, ~, ~]  = NumFacsDFM_FAVAR(X, MaxFac);
    saveas(gcf,['../Figures/PostVolckNumFacs.pdf']);
    K = nfact(:,2500); %% 2500'th entry seems stable, K = number of factors;
    [NumLags, IC_table] = DetermineLagLengthFAVAR(K, faclag, MaxLag, Y, X_st, slowcode, det,type);
    plag = NumLags(2);         % plag is number of lags in the VAR part
else
    K = 3;
    plag = 3;
end
clearvars vnames
for i = 1:K
   vnames{i} = strcat("Factor ", char(string(i)));
end
vnames = [vnames';namesY];
%% FAVAR
if K == 3 && plag == 3
    if Cholesky == 1
%-----------------------Reduced form estimation----------------------------
[FAVAR_Chol, FAVARopt_Chol] = EstimateReducedFAVAR(X_st, Y, faclag, plag, K, det, slowcode, type, namesXY,MLGARCH,1);
FAVAR_Cusum(FAVAR_Chol,DATES,namesXY,K,['PostVolck_Cusum',firstdstr,num2str(MLGARCH)],1);
[FAVARPPCA, ~] = EstimateReducedFAVAR(X_st, Y, 1, plag, K, det, slowcode, 1, namesXY,MLGARCH,Cholesky);
RMSEPost = sqrt(mean(FAVAR_Chol.residuals.^2));
if MLGARCH == 1
    FAVARopt_Chol.vnames = vnames;
PlotFilteredVolatilityFAVAR(FAVAR_Chol, FAVARopt_Chol, DATES, ['Filtered_Vola_subsample',firstdstr]);
PlotResidualFAVARGarch(FAVAR_Chol,FAVARopt_Chol,DATES,'FAVAR_post_Residuals');
end
qqplot_FAVARResid(FAVAR_Chol, FAVARopt_Chol, ['QQ_19832019',firstdstr,num2str(MLGARCH)]);
PlotEigs(FAVAR_Chol,['EigenVals19832019',firstdstr]);
[GCtest]  = FAVARGrangerCausalityTest(FAVAR_Chol.ENDO, 1, vnames);
disp(round(GCtest.Variables,3))
Plot_FAVARfilter(FAVAR_Chol, FAVARopt_Chol, DATES, ['Filtered_FAVAR_Subsample',firstdstr]);
%table2latex(GCtest,'../Tables/GCTest.tex');
%% R2-gain over PCA
figure
set(0,'DefaultAxesColorOrder',[0 0 0]);
DFMR2s = FAVAR_Chol.DFM.R2s.Variables;
PCAR2s = FAVARPPCA.DFM.R2s.Variables;
histogram(DFMR2s(1:end-2) - PCAR2s(1:end-2),'BinWidth',0.025,'Normalization','probability')
hold on
mu = xline(mean(DFMR2s(1:end-2) - PCAR2s(1:end-2)),'--','LineWidth',1.5);
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
legend(mu, ['Mean = ', num2str(round(mean(DFMR2s(1:end-2) - PCAR2s(1:end-2)),3))],'Interpreter','latex','FontSize',14)
ylabel('Relative frequency'); xlabel('Marginal Gain: $R^2_{DFM} - R^2_{PCA}$')
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/HistogramR2PostVolck.pdf']);
%% Total R^2 of the panel histofram
figure
set(0,'DefaultAxesColorOrder',[0 0 0]);
DFMR2s = FAVAR_Chol.DFM.R2s.Variables;
histogram(DFMR2s(1:end-2),'BinWidth',0.025,'Normalization','probability')
hold on
mu = xline(mean(DFMR2s(1:end-2)),'--','LineWidth',1.5);
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
legend(mu, ['Mean = ', num2str(mean(DFMR2s(1:end-2)))],'Interpreter','latex','FontSize',14)
ylabel('Relative frequency'); xlabel('Marginal Gain: $R^2_{DFM} - R^2_{PCA}$')
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/HistogramR2TotalPostVolck.pdf']);
%%
R2tab = table(PCAR2s(var_numbers), DFMR2s(var_numbers));
R2tab.Properties.RowNames = namesXY(var_numbers); R2tab.Properties.VariableNames = ["PCA", "DFM"]
diff = sum([DFMR2s(var_numbers), -PCAR2s(var_numbers)],2);
    end
end
%%
if Cholesky == 1
[FAVAR, FAVARopt] = EstimateReducedFAVAR(X_st, Y, faclag, plag, K, det, slowcode, type, namesXY,MLGARCH,Cholesky);
FAVARopt.vnames = vnames;
Plot_FAVARfilter(FAVAR, FAVARopt, DATES, ['Filtered_FAVAR_Subsample',firstdstr]);
R2s = FAVAR.DFM.R2s(var_numbers,1);
R2s
%% FAVAR - Cholesky Short-Run Restrictions:
%---------------------------Cholesky SVAR----------------------------------
% Set options some options for IRF calculation and bands
FAVARCholopt = GenerateOptionsFAVAR(1, ConfLevel, nsteps, vnames);
%-----------------------Impulse Response Functions-------------------------
[IRF_Chol, FAVAR_Chol] = VARir(FAVAR,FAVARCholopt);
% Compute error bands
[IRFINF_Chol,IRFSUP_Chol,IRFMED_Chol,IRFINF_Cholx,IRFSUP_Cholx,IRFMED_Cholx] = FAVARirband(FAVAR_Chol,FAVARCholopt);
save("IRFMedCholPost", "IRFMED_Chol")
%---------------------------------VAR IRF----------------------------------
nlagstr = num2str(FAVAR_Chol.nlag);
name  = ['Chol',nlagstr,'Lags19832019',firstdstr];
FAVAR_Resp(FAVAR_Chol, FAVARCholopt, IRFMED_Chol, IRFINF_Chol, IRFSUP_Chol,name, FFR1d,1);
%---------------------------------DFM IRF----------------------------------
prefix = ['Chol',nlagstr,firstdstr];
IRF_DFM(FAVAR, FAVARCholopt, IRFMED_Cholx, IRFINF_Cholx, IRFSUP_Cholx,  tcodes, var_numbers, var_names, prefix, Order, FFR1d)
%----------------------Historical Decomposition----------------------------
if HD == 1
Modelname = ['3lagChol1983',firstdstr];
FAVARCholopt.pick = 1;
set(0,'DefaultAxesColorOrder','remove');
HD = VARhd(FAVAR_Chol);
VARhdplot(HD,FAVARCholopt, FAVAR_Chol, DATES, Modelname);
end
var_names  = [var_names;namesY];
%--------------------------------DFM FEVD----------------------------------
set(0,'DefaultAxesColorOrder','remove');
FAVARCholopt.nsteps = 24;
[FEVD, VAR, FAVARFEVD] = FAVARfevd(FAVAR,FAVARCholopt);
figure
for i = 1:size(var_numbers, 2)
index = var_numbers(i);
    subplot(4,3,i)
Plootly = FAVARFEVD(:,:,index);
area(Plootly);
title(var_names{i});
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
end
set(gcf, 'unit', 'inches');
figure_size =  get(gcf, 'position');
lg = legend(FAVARCholopt.vnames,'orientation','horizontal','location','southoutside','interpreter','latex');
set(lg,'Units','inches',...
    'Position',[3.5 0.1 1.0 0.23],...
    'Interpreter','latex',...
    'FontSize',14);
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/FEVD',prefix,nlagstr,'.pdf']);
else
%% Blanchard & Quah Long Run Restrictions
prefix = 'bq';
Y = ydata_st;
%---------------------------Long Run SVAR----------------------------------
[FAVAR, ~] = EstimateReducedFAVAR(X_st, Y, faclag, plag, K, 0, slowcode, type,namesXY,MLGARCH,Cholesky);
FAVARbqopt = GenerateOptionsFAVAR(2, 68, nsteps, vnames);
FAVARbqopt.standardized = 1;
FAVARbqopt.snames = [vnames(1);"Aggregate Supply"; "Spending";vnames(2:K);"Expansionary Monetary Policy"]; 
FAVARbqopt.vnames = [vnames(1);"IP"; "FFR";vnames(2:K);"CPI"]; 
FAVARbqopt.NegativeShockIndex = [];
plotCOVFAVAR; % Seperate script plotting the conditional covariance
% Compute IRF
[IRF_bq, FAVAR_bq] = VARir(FAVAR,FAVARbqopt);
% Compute error bands
[IRFINF_bq,IRFSUP_bq,IRFMED_bq,INFx_bq,SUPx_bq,MEDx_bq] = FAVARirband(FAVAR_bq,FAVARbqopt);
if K == 3 || plag == 3
    save("IRFMedbqPost", "IRFMED_bq")
end
%---------------------------Long Run IRF-----------------------------------
nlagstr = num2str(FAVAR_bq.nlag);
name  = [prefix,nlagstr,'Lags19832019',firstdstr,'K',num2str(K)];
FAVAR_Resp(FAVAR_bq, FAVARbqopt, IRFMED_bq, IRFINF_bq, IRFSUP_bq,name, FFR1d,false);
%-------------------------Long Run IRF DFM---------------------------------
IRF_DFM(FAVAR_bq, FAVARbqopt, MEDx_bq, INFx_bq, SUPx_bq, tcodes, var_numbers, var_names, prefix, Order, FFR1d)
var_names = [var_names;namesY];
%%
if plag == 3 && K == 3
%---------------------------Long Run FEVD----------------------------------
%%
clearvars INF SUP MED INFx SUPx MEDx
FAVARbqopt.nsteps = 24;
tic
[FEVD, VAR, FAVARFEVD] = FAVARfevd(FAVAR_bq,FAVARbqopt);
[INF,SUP,MED,INFx,SUPx,MEDx] = FAVARfevdband(FAVAR_bq,FAVARbqopt);
toc
set(0,'DefaultAxesColorOrder','remove');
%%
figure
for i = 1:size(var_numbers, 2)
index = var_numbers(i);
    subplot(4,3,i)
Plootly = MEDx(:,:,index);
area(Plootly);
ylim([0,1]);
title(var_names(i));
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
end
% set unit for figure size to inches
set(gcf, 'unit', 'inches');
lg = legend(FAVARbqopt.snames,'orientation','horizontal','location','southoutside');
set(lg,'Units','inches',...
    'Position',[3.5 0.1 1.0 0.23],...
    'Interpreter','latex',...
    'FontSize',14);
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/FEVD',prefix,num2str(plag),'Lags.pdf']);
%%
% Generate tables:
for i=1:size(INFx,1)
    for j=1:size(INFx,2)
        for k=1:size(INFx,3)
            CIband{i,j,k} = string(['(',num2str(round(INFx(i,j,k)*100,1)),', ',num2str(round(SUPx(i,j,k)*100,1)),')']);
        end
    end
end
%%
horizons = [1,2,3,6,12,24];
CIband = CIband(horizons,:,var_numbers);
MEDFEVD = round(MEDx(horizons,:,var_numbers),4);
tablenames = ["IPFEVD.tex", "FFRFEVD.tex", "CPIFEVD.tex"];
caption = ["FEVD: Industrial Production", "FEVD: Federal Funds Rate","FEVD: Consumer Price Index"];
for k = 1:3
ncol = size(CIband,2);
nrow = size(CIband,1);
col_spec = ['l'];
    for c = 1:ncol
        col_spec = [col_spec 'c']; 
    end
col_names    = strjoin(["\textit{horizons}";FAVARbqopt.snames], ' & ');
row_names    = strsplit(string(num2str(horizons)));
filename = strjoin(['../Tables/',tablenames(k)],'');
fileID = fopen(filename, 'w');
fprintf(fileID, '\\begin{table}');
fprintf(fileID, strjoin(['\\caption{', caption(k), '}']));
fprintf(fileID,'\\begin{adjustbox}{max width=\\linewidth} \n');
fprintf(fileID, '\\begin{tabular}{%s}\n', col_spec);
fprintf(fileID, '\\toprule \n');
fprintf(fileID,' & \\multicolumn{2}{c}{\\textit{Supply components}} & \\multicolumn{4}{c}{\\textit{Demand Components}} \\\\ \n');
fprintf(fileID,' \\cmidrule(l){2-3} \\cmidrule(l){4-7} \n');
fprintf(fileID, '%s \\\\ \n', col_names);
fprintf(fileID, '\\midrule \n');
for i = 1:nrow
   temp{1,ncol}  = []; 
   temp1{1,ncol} = []; 
   for col = 1:ncol
      value  = MEDFEVD(i,col,end-3+k);
      value1 = CIband{i,col,end-3+k};
      temp{1,col}  = string(num2str(value*100));
      temp1{1,col} = value1;
   end
    fprintf(fileID, '%s \\\\ \n', strjoin([row_names(i) string(temp)], ' & '));
    fprintf(fileID, '%s \\\\ \n', strjoin([" ", string(temp1)], ' & '));
    if i < nrow
    fprintf(fileID,'\\addlinespace \n');
    end
end
 fprintf(fileID, '\\bottomrule \n');
    fprintf(fileID, '\\end{tabular} \n');
    fprintf(fileID, '\\end{adjustbox}\n');
     fprintf(fileID, '\\begin{minipage}{\\textwidth} \n');
     fprintf(fileID, '{\\footnotesize\n');
     fprintf(fileID, '\\begin{enumerate}\n');
     fprintf(fileID, '\\vspace{2mm}\n');
     fprintf(fileID,'    \\item[1] Forecast error variance decomposition of the FAVAR; brackets below the point estimate shows one standard deviation error bands\n');
     fprintf(fileID,'\\end{enumerate}}\n');
     fprintf(fileID,'     \\end{minipage}\n');
     fprintf(fileID, '\\end{table}');
    fclose(fileID);
end
end
end % End of structural analysis
%% OOS, Forecast plot
if K == 3 && plag == 3;
set(0,'DefaultAxesColorOrder','default');
clearvars Fitntest Datetesttest

figure;
h = 12;
WinSize = 120;
start = 200;
tests = [start:25:400, 430-h-1];
iter = 0;
    for i = tests
        iter = iter + 1;
        [~, Fitntest(:,:,iter), ~, Datetesttest(iter,:)] = RollingWindowFAVARComp(X_st, Y,faclag, plag,...
            K, 0, slowcode,DATES,WinSize,i,type,h,true,namesXY);
    end
for Ser = 1:3
subplot(3,1,Ser)
obsplot = plot(DATES(start-1:430), Y(start-1:430,Ser),'-','LineWidth',1,'Color',[.5 .5 .5]);
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
title(namesY{Ser});
for i = 1:size(tests,2)
    hold on
    if i == 1 | Ser == 1
        fcplot = plot(Datetesttest(i,:), Fitntest(:,Ser,i),'k-','LineWidth',1.25);
    else
    plot(Datetesttest(i,:), Fitntest(:,Ser,i),'k-','LineWidth',1.25)
    end
end
end
% set unit for figure size to inches
set(gcf, 'unit', 'inches');
% get the original size of figure before the legends are added
figure_size =  get(gcf, 'position');
lg = legend([obsplot, fcplot], {'Observation', 'Forecast'},'orientation','horizontal','location','southoutside','interpreter','latex');
set(lg,'Units','inches',...
    'Position',[3.5 0.1 1.0 0.23],...
    'Interpreter','latex',...
    'FontSize',14);
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/OoSForecasts1983.pdf']);
end