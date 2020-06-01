clear 
addpath(genpath('Data'));
addpath(genpath('Utils'));
clear all; clc;
set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
Order = [1, 2, 3];
FFR1d = 1;
PrepdataFull
diffFFR = [NaN; ydata(2:end,3) - ydata(1:end-1,3)];
%% Util script fetching 2 of our 3 series
url = 'https://fred.stlouisfed.org/';
c = fred(url);
%% Retrieve seasonally adjusted series
series = [{'CPIAUCSL', 'INDPRO'}];
startdate = '01/01/1959';
enddate = '11/01/2019';
Y = [];
for i=1:length(series)
    d = fetch(c,series(i),startdate,enddate);
    dateY = d.Data(:,1);
    Y = [Y, d.Data(:,2)];
end
%% 
logD_Y = [NaN, NaN; log(Y(2:end,:))-log(Y(1:end-1,:))];
%%
logD_Y = [logD_Y(4:end,:), diffFFR];
datelogY = dateY(2:end,:);
Y = Y(1:end-3,:);
Y = [Y, ydata(:,3)];
%%
for i = 1:3
    [~,pValue(i),stat(i),cValue(i),~] = adftest(Y(:,i),'Lags',3);
    [~,logpValue(i),logstat(i),logcValue(i),~] = adftest(logD_Y(:,i),'Lags',3);
end
tableADF = table([stat;cValue;logstat;logcValue]);
%%
tableADF.Properties.RowNames = ["StatLevel", "CritLevel", "StatLog", "CritLog"];
tableADF = splitvars(tableADF);
tableADF.Properties.VariableNames = ["CPI", "IP", "FFR"]
% %table2latex(tableADF,'../Tables/ADF.tex');
%%
set(0,'DefaultAxesColorOrder',rgb([0, 0, 0]));
subplot(2,1,1)
plot(datetime(dateY(4:end),'ConvertFrom','datenum'), Y(:,1));
title('CPI','Interpreter','latex','FontSize',16)
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
ylabel('$CPI$','Interpreter','latex')
ylim([min(Y(:,1)) - std(diff(Y(:,1))), max(Y(:,1))+std(diff(Y(:,1)))]);
recessionplot_RMJ("k",0.1); 
subplot(2,1,2)
plot(datetime(datelogY(3:end),'ConvertFrom','datenum'), logD_Y(:,1))
title('Inflation Rate','Interpreter','latex','FontSize',16)
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
ylabel('$\Delta \ln \left( CPI \right)$','Interpreter','latex')
xlabel('Date')
ylim([min(logD_Y(:,1)) - std(logD_Y(:,1)), max(logD_Y(:,1))+std(logD_Y(:,1))]);
recessionplot_RMJ("k",0.1); 

set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/InflationSeries.pdf']);
%%
set(0,'DefaultAxesColorOrder',rgb([0, 0, 0]));
subplot(2,1,1)
plot(datetime(dateY(4:end),'ConvertFrom','datenum'), Y(:,2));
title('Industrial Production','Interpreter','latex','FontSize',16)
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
ylabel('$IP$','Interpreter','latex')
ylim([min(Y(:,2)) - std(diff(Y(:,2))), max(Y(:,2))+std(diff(Y(:,2)))]);
recessionplot_RMJ("k",0.1); 
subplot(2,1,2)
plot(datetime(datelogY(3:end),'ConvertFrom','datenum'), logD_Y(:,2))
title('Growth rate of the industrial production','Interpreter','latex','FontSize',16)
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
ylabel('$\Delta \ln \left(IP\right)$','Interpreter','latex')
xlabel('Date')
ylim([min(logD_Y(:,2)) - std(logD_Y(:,2)), max(logD_Y(:,2))+std(logD_Y(:,2))]);
recessionplot_RMJ("k",0.1); 

set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/IndustrialProductionSeries.pdf']);
%%
writematrix(logD_Y, "seasonallyAdjData.xlsx")
%%
FFR1d = 1;
PrepdataFull
%%
Q = readmatrix('FEDFUNDS-3.xls');
M = readmatrix('FEDFUNDS-2.xls');
%% 
Mdate = datetime(M(:,1),'ConvertFrom','datenum');
Mdate = Mdate + 365.25 * 1900;
Qdate = datetime(Q(:,1),'ConvertFrom','datenum');
Qdate = Qdate + 365.25 * 1900;
ELB = [datenum(2009,01,05), datenum(2016,01,04)];
%%
SampleInd = Mdate >= DATES(1) & Mdate <= DATES(end);
SampleIndQ = Qdate >= DATES(1) & Qdate <= DATES(end);
%%
figure;
subplot(3,1,1)
plot(Mdate(SampleInd), M(SampleInd,2),'LineWidth',1.25);
title('Effective Federal Funds Rate, Monthly','FontSize',16,'Interpreter','latex');
recessionplot_RMJ("k",0.075,"recession",ELB); 
xline(datetime([1984, 1, 1]),':', 'LineWidth',1.25)
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
%recessionplot_RMJ("k",0.1);
subplot(3,1,2)
plot(Qdate(SampleIndQ), Q(SampleIndQ,2),'LineWidth',1.25);
title('Effective Federal Funds Rate, Quarterly','FontSize',16,'Interpreter','latex');
recessionplot_RMJ("k",0.075,"recession",ELB); 
xline(datetime([1984, 1, 1]),':', 'LineWidth',1.25)
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
subplot(3,1,3)
plot(DATES, ydata(:,3),'LineWidth',1.25);
title('Policy Rate, Monthly','FontSize',16,'Interpreter','latex');
ylim([min(ydata(:,3))-0.5*std(ydata(:,3)), max(ydata(:,3))+0.5*std(ydata(:,3))]);
yline(0,'--','LineWidth',1.1);
xline(datetime([1984, 1, 1]),':', 'LineWidth',1.25)
recessionplot_RMJ("k",0.075,"recession",ELB); 
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
%recessionplot_RMJ("k",0.5); 
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
%%saveas(gcf,['../Figures/FFRMQL.pdf']);
%%
FFR1diff = ydata(2:end,3) - ydata(1:end-1,3);
set(0,'DefaultAxesColorOrder',rgb([0, 0, 0]));
plot(DATES(2:end), FFR1diff)
recessionplot_RMJ("k",0.075,"recession",ELB); 
recessionplot_RMJ("k",0.1); 
xline(datetime([1984, 1, 1]),':', 'LineWidth',1.25)
title('$\Delta$Policy Rate, Monthly','FontSize',16,'Interpreter','latex');
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
xline(datetime([1984, 1, 1]),':', 'LineWidth',1.25)
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
%saveas(gcf,['../Figures/FFR1d.pdf']);
%%
iter = 0;
for i=100:500
    iter = iter+1;
[~,~,statInf,cValueInf,~] = kpsstest(ydata_st(1:i,3),'lags',1:12);
[~,~,statSup,cValuesup,~] = kpsstest(ydata_st(i+1:end,3),'lags',1:12);
breakpoint(iter,:) = i;
statInfmat(iter,:) = statInf; statSupmat(iter,:) = statSup;
cValueInfmat(iter,:) = cValueInf;cValuesupmat(iter,:) = cValuesup;
end
iter = 0;
for i=100:500
    iter = iter+1;
[~,~,statInf,cValueInf,~] = adftest(ydata_st(1:i,3),'lags',1:12);
[~,~,statSup,cValuesup,~] = adftest(ydata_st(i+1:end,3),'lags',1:12);
breakpoint_ADF(iter,:) = i;
statInfmat_ADF(iter,:) = statInf; statSupmat_ADF(iter,:) = statSup;
cValueInfmat_ADF(iter,:) = cValueInf;cValuesupmat_ADF(iter,:) = cValuesup;
end
%%
set(0,'DefaultAxesColorOrder',rgb([0, 0, 0]));
subplot(2,1,1)
plot(DATES(100:500),statInfmat_ADF(:,12),':','LineWidth',1.5)
hold on
plot(DATES(100:500),statInfmat_ADF(:,6),'-.','LineWidth',1.5) 
hold on
plot(DATES(100:500),statInfmat_ADF(:,1),'--','LineWidth',1.5) 
hold on
plot(DATES(100:500), cValueInfmat_ADF(:,12),'-','LineWidth',1.25)
ylabel('ADF-statistic','interpreter','latex');
xlabel('Ultimo date','interpreter','latex');
yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
subplot(2,1,2)
plot(DATES(end-400:end),statSupmat_ADF(:,12),':','LineWidth',1.25) %% plotter over alle iterationerne med 1 lag
hold on
plot(DATES(end-400:end),statSupmat_ADF(:,6),'-.','LineWidth',1.25) 
hold on
plot(DATES(end-400:end),statSupmat_ADF(:,1),'--','LineWidth',1.25) 
hold on
plot(DATES(end-400:end), cValuesupmat_ADF(:,12),'-','LineWidth',1.1)
legend('$p=12$', '$p=6$', '$p=1$','Critical value, $\alpha=0.05$','interpreter','latex','Location','best','FontSize',30);
ylabel('ADF-statistic','interpreter','latex');
xlabel('Primo date','interpreter','latex');
yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
%legend('$p=12$', '$p=6$', '$p=1$','Critical value, $\alpha=0.05$','interpreter','latex','Location','best','FontSize',14);
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/IterADF1d.pdf']);
%%
set(0,'DefaultAxesColorOrder',rgb([0, 0, 0]));
subplot(2,1,1)
plot(DATES(100:500),statInfmat(:,12),':','LineWidth',1.5) %% plotter over alle iterationerne med 1 lag
hold on
plot(DATES(100:500),statInfmat(:,6),'-.','LineWidth',1.5) 
hold on
plot(DATES(100:500),statInfmat(:,1),'--','LineWidth',1.5) 
hold on
plot(DATES(100:500), cValueInfmat(:,12),'-','LineWidth',1.25)
ylabel('KPSS-statistic','interpreter','latex');
xlabel('Ultimo date','interpreter','latex');
yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
legend('$p=12$', '$p=6$', '$p=1$','Critical value, $\alpha=0.05$','interpreter','latex','Location','best','FontSize',30);
subplot(2,1,2)
plot(DATES(end-400:end),statSupmat(:,12),':','LineWidth',1.25) %% plotter over alle iterationerne med 1 lag
hold on
plot(DATES(end-400:end),statSupmat(:,6),'-.','LineWidth',1.25) 
hold on
plot(DATES(end-400:end),statSupmat(:,1),'--','LineWidth',1.25) 
hold on
plot(DATES(end-400:end), cValuesupmat(:,12),'-','LineWidth',1.1)
ylabel('KPSS-statistic','interpreter','latex');
xlabel('Primo date','interpreter','latex');
yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
%legend('$p=12$', '$p=6$', '$p=1$','Critical value, $\alpha=0.05$','interpreter','latex','Location','best','FontSize',14);
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
% saveas(gcf,['../Figures/IterKPSS1d.pdf']);
%% Retrieve seasonally adjusted series
url = 'https://fred.stlouisfed.org/';
c = fred(url);
series = [{'INDPRO', 'CPIAUCSL'}];
startdate = '12/01/1958';
enddate = '12/01/2019';
Y = [];
for i=1:length(series)
d = fetch(c,series(i),startdate,enddate);
dateY = d.Data(:,1);
Y = [Y, d.Data(:,2)];
end
%%
FFR1d = 0;
PrepdataFull 
New = [Y(4:end,:), ydata(:,end)];
%% 
[h,pValue,stat,cValue,mles] = jcitest(New(1:300,:),'model','H1','lags',3);
CE = 0:2;
stat = table2array(stat)';
cVal = table2array(cValue)';
pValue = table2array(pValue)';
TracePre = table(CE', stat, cVal,pValue);
%%
[h,pValue,stat,cValue,mles] = jcitest(New(301:end,:),'model','H1','lags',3);
CE = 0:2;
stat = table2array(stat)';
cVal = table2array(cValue)';
pValue = table2array(pValue)';
TracePost = table(CE', stat, cVal,pValue);
%%
[h,pValue,stat,cValue,mles] = jcitest(New,'model','H1','lags',3);
CE = 0:2;
stat = table2array(stat)';
cVal = table2array(cValue)';
pValue = table2array(pValue)';
TraceFull = table(CE', stat, cVal,pValue);
%%
table2latex(TraceFull,'FullJohansen.tex')
table2latex(TracePost,'PostJohansen.tex')
table2latex(TracePre,'PreJohansen.tex')
%% Level Summary Stats:
FullSum = GenerateSummaryStats(New, namesY);
PreSum = GenerateSummaryStats(New(1:300,:), namesY);
PostSum = GenerateSummaryStats(New(301:end,:), namesY);
%%
PrepdataFull
%% Level Summary Stats:
FullSum = GenerateSummaryStats(ydata, namesY);
PreSum  = GenerateSummaryStats(ydata(1:300,:), namesY);
PostSum = GenerateSummaryStats(ydata(301:end,:), namesY);
%%
table2latex(FullSum,'FullSum.tex');
table2latex(PreSum,'PostSum.tex');
table2latex(PostSum,'PreSum.tex');
%% 
for i=1:12 % 
[h,pValue,stat,cValue,reg] = adftest(ydata_st(:,3),'lags',i);
BIC(i) = reg.BIC;
end
[~, OptLag] = min(BIC);
[h,pValue,stat,cValue,reg] = adftest(ydata(1:289,3),'lags',OptLag);
% autocorr_RMJ(reg.res);
% pValue
%PreSum  = GenerateSummaryStats(ydata(1:288,:), namesXY(129:131));
%PostSum = GenerateSummaryStats(ydata(289:end,:), namesXY(129:131));
%%
clear all; clc; close all;
% Settings
signRest         = 0;         % Aux switch indicating whether or not sign-restrictions should be applied
FFR1d            = 1;         % Switch to determine whether FFR should be assumed I(1) or levels I(0)
Pre2009          = 0;         % Aux switch, 0 = 1959:2019, 1 = 1959:2008
type             = 2;         % Switch 1 = PCA extraction of factors, 2 = Kalman Smoothing of factors. NB. Kalman is computationally very heavy.
nsteps           = 24;        % Number of steps IRF
ConfLevel        = 68;        % Confidence level percentage for bootstrap
Information_Crit = 0;         % Switch 1: Lag length and numfacs by IC
Order            = [1,2,3];      % 1: Inflation, 2: Industrial Prod, 3: FFR
                                 % Cholesky order [1,2,3];
                                 % BQ order: [2,3,1];
ELB = [datenum(2009,01,05), datenum(2016,01,04)];
var_numbers = [1, 19, 24, 66, 73, 80, 128, 129, 130 131];
PrepData;
%%
set(0,'DefaultAxesColorOrder',rgb([0, 0, 0]));
figure;
plot(DATES, ydata(:,3));
%title('$\Delta$Policy Rate','Interpreter','latex','FontSize',16);
recessionplot_RMJ("k",0.075,"recession",ELB); 
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/FFR1dS.pdf']);
