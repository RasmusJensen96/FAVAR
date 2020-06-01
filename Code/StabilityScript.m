%-------------------------------------------------------------------------%
%                   Structural Break tests and plots                      %
%-------------------------------------------------------------------------%
% This script does data preparation and conducts structural break tests for
% the FAVAR specification as identified in main. It uses the forecast
% equation for each observable variable from the reduced form VAR and runs
% a Chow-test and two Cusum-tests, for each of the sample specifications:
% 1959-2019, and 1983-2019.
%
% Rasmus M. Jensen 2020.
clear all; clc; close all;
addpath(genpath('Data'));
addpath(genpath('Utils'));

rng(260196,'twister');
set(0,'DefaultAxesColorOrder',[0 0 0]);
FFR1d = 1;
Order = [1,2,3];
faclag = 3;
PrepdataFull;
PV = [datenum(1979,08,06), datenum(1987,08,11)]; % Paul Volcker office period.
figure;
for i = 1:3
subplot(3,1,i)
plot(DATES, ydata_st(:,i))
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
ylim([min(ydata_st(:,i))- 0.25*std(ydata_st(:,i)), max(ydata_st(:,i)) + 0.25* std(ydata_st(:,i))]);
if i == 3
    xlabel('Date','interpreter','latex','FontSize',14)
end
ylabel('Standardized value','interpreter','latex','FontSize',14)
title(namesXY{125+i,:},'interpreter','latex','FontSize',16);
recessionplot_RMJ("k",0.20);
recessionplot_RMJ("k",0.075,"recession",PV); 
end 
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/StandardizeSeriesFullSample.pdf']);
% Chow-test:
t = 22;
J = 700;
[FAVAR, FAVARopt] = EstimateReducedFAVAR(x2_st, ydata_st, 3, 3, 3, 0, slowcode, 2, namesXY,1,1,0,0);
FY_GARCH = FAVAR.fit + FAVAR.residuals;
FY_Gaussian = FAVAR.fit + FAVAR.residualsNonFilter;
[ChowStat, cValue] = IterChow(FY_GARCH,3,3, t, J, DATESnm,1,1,'Chow_Full_sampleFAVAR_GARCH');
[ChowStat, cValue] = IterChow(FY_Gaussian,3,3, t, J, DATESnm,1,1,'Chow_Full_sampleFAVAR');
%% Cusum-test
name = ["Gaussian";"GARCH"];
for jj = 1:2
eval(['FY = FY_',char(name(jj))]);
plag = 3;
K    = 3;
N = size(FY,2);
RHS = FY(plag+1:end,K+1:N); % Observables
%RHS = FY(plag+1:end,1:3); % Factors
LHS = FY(plag:end-1,:);
for L = 1:plag-1
    LHS = [LHS, FY(plag-L:end-L-1,:)];
end
%% 
Testtype = [{'cusum'},{'cusumsq'}];
for j=1:2
        for i=1:3
            [~,~,Stat,~,~ ,c1,c2] = cusumtest_modded(LHS, RHS(:,i),'Test',Testtype(j));
            Cusum_Stat(:,i,j) = Stat;
            Cusum_upper(:,i,j) = c1;
            Cusum_lower(:,i,j) = c2;
end
end
%%
index= reshape(1:6, 2, 3)';
Testtype = ["residuals", "squared residuals"];
    figure       
    t = 25;
iter = 0;
for j = 1:2
    for i=1:3
        iter = iter+1;
        subplot(3,2,index(iter));   
        plot(DATES(t:end-1,1), Cusum_upper(:,i,j),'LineStyle','--','LineWidth',1.25);
        hold on
        plot(DATES(t:end,1), Cusum_Stat(:,i,j),'LineWidth',1.5);
        hold on 
        plot(DATES(t:end-1,1), Cusum_lower(:,i,j),'LineStyle','--','LineWidth',1.25);
        %title(namesXY(size(namesX,1)+i),'FontSize',16,'interpreter','latex'); 
        title(namesXY(size(x2,2)+i),'FontSize',16,'interpreter','latex'); 
       if i == 1
           title({['Cumulative sum of ' Testtype{j}], 'Inflation'},'FontSize',16,'interpreter','latex');
          legend('Critical Value, $5\%$','Test statistic ','interpreter','latex','FontSize',14,'location','southeast'); 
       else
           title(namesXY(size(x2,2)+i),'FontSize',16,'interpreter','latex'); 
       end
       betweenShading = [DATES(t:end-1)', fliplr(DATES(t:end-1)')];
       CurvatureShade = [Cusum_lower(:,i,j)', fliplr(Cusum_upper(:,i,j)')];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
       xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
       if j == 2
       ylim([-0.1268 1]); % Recessionband correction when auto-exporting
       elseif j==1
           ylim([-100 120]);
       end
       %recessionplot;
       recessionplot_RMJ("k",0.20);
       recessionplot_RMJ("k",0.075,"recession",PV); 
    end
end
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/CusumFactor_Full_sample',char(name(jj)),'.pdf']);
end
%%
for jj = 1:2
    eval(['FY = FY_',char(name(jj))]);
    N = size(FY,2);
    NF = 3;
    NL = plag;
    NY = N - NF;
    RHS = FY(NL+1:end,1:3);
    LHS = FY(NL:end-1,:);
    for L = 1:plag-1
        LHS = [LHS, FY(NL-L:end-L-1,:)];
    end
%% 
Testtype = [{'cusum'},{'cusumsq'}];
for j=1:2
        for i=1:3
[~,~,Stat,~,~ ,c1,c2] = cusumtest_modded(LHS, RHS(:,i),'Test',Testtype(j));
Cusum_Stat(:,i,j) = Stat;
Cusum_upper(:,i,j) = c1;
Cusum_lower(:,i,j) = c2;
end
end
%%
namesF = ["Factor 1", "Factor 2", "Factor 3"];
index= reshape(1:6, 2, 3).';
Testtype = ["residuals", "squared residuals"];
    figure       
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    namesXY(end,:) = {'Policy Rate'};
    t = 25;
iter = 0;
for j = 1:2
    for i=1:3
        iter = iter+1;
        subplot(3,2,index(iter));   
        plot(DATES(t:end-1,1), Cusum_upper(:,i,j),'LineStyle','--','LineWidth',1.25);
        hold on
        plot(DATES(t:end,1), Cusum_Stat(:,i,j),'LineWidth',1.5);
        hold on 
        plot(DATES(t:end-1,1), Cusum_lower(:,i,j),'LineStyle','--','LineWidth',1.25);
       if j == 1 
           mini = min(min(Cusum_Stat(:,i,j)),min(Cusum_lower(:,i,j)))-20;
           maxi = max(max(Cusum_Stat(:,i,j)),max(Cusum_upper(:,i,j)))+20;
           ylim([mini maxi]);
       else
           mini = min(min(Cusum_Stat(:,i,j)),min(Cusum_lower(:,i,j)))-0.20;
           maxi = max(max(Cusum_Stat(:,i,j)),max(Cusum_upper(:,i,j)))+0.20;
           ylim([mini maxi]);
       end
        title(namesF(i),'FontSize',16,'interpreter','latex'); 
       if i == 1
           title({['Cumulative sum of ' Testtype{j}], 'Factor 1'},'FontSize',16,'interpreter','latex');
          legend('Critical Value, $5\%$','Test statistic ','interpreter','latex','FontSize',14,'location','Southeast'); 
       else
           title(namesF(i),'FontSize',16,'interpreter','latex'); 
       end
       yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
        xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
       betweenShading = [DATES(t:end-1)', fliplr(DATES(t:end-1)')];
       CurvatureShade = [Cusum_lower(:,i,j)', fliplr(Cusum_upper(:,i,j)')];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       recessionplot_RMJ("k",0.20);
       recessionplot_RMJ("k",0.075,"recession",PV);
    end
end
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/CusumFactors_Full_sample',char(name(jj)),'.pdf']);
end
 %%
%-----------------------------Factor Equation------------------------------
clearvars -except Order PV
FFR1d = 1;
PrepData;
K = 3;               % Number of Factors
plag = 3;            % plag is number of lags in the VAR part
X = x2_st;
X_st = x2_st;
Y = ydata_st;
%%
%-----------------------------Factor Equation------------------------------
% % Extract principal components from (standardized) X
% X_st = x2_st; 
% [F0,~]=extract(X_st,K); % F0 are the factors, Lf are the loadings
% %Now rotate the factor space as in Bernanke, Boivin and Eliasz (2005)
% slowindex = find(slowcode==1)';
% xslow = X(:,slowindex);
% [Fslow0,~] = extract(xslow,K);
% Fr0 = facrot(F0,Y(:,end),Fslow0);
F0 = DFM_extraction(X_st, 2, K, 1);
% Rotating factor space as in Bernanke, Boivin and Eliasz (2005)
slowindex = find(slowcode==1)';
xslow = X_st(:,slowindex);
%Fslow0 = DFM_extraction(xslow, 1, K, 1);
% [Fslow0,~] = extract(xslow,K);
%Fr0 = facrot(F0,Y(:,end),Fslow0);
Fr0 = F0; % (Fr0-mean(Fr0))./std(Fr0)
FY=[Fr0,Y];  % the extracted factors, plus Y (infl,unemp and interest)
vnames = ["Factor 1", "Factor 2", "Factor 3"];
%% Cond. Mean model Factors: ACF, PACF
%---------------------Cond mean model for factors--------------------------
F = Fr0;
index = reshape(1:6,2,3).';
iter = 0;
figure       
for i=1:numel(index)
    iter = iter + 1;
    subplot(3,2,index(iter))
    if i == 1 | i==2 | i==3;
        autocorr_RMJ(F(:,i));
        title(vnames(i),'Interpreter','latex','FontSize',14)
        yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
        xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
        if i == 3
            xlabel('Lags');
        end
    else
        parcorr_RMJ(F(:,i-3));
        title(vnames(i-3),'Interpreter','latex','FontSize',14)
        yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
        xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
        if i == 6
            xlabel('Lags');
        end
    end
end
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/AutoCorrFactor_PCA.pdf']);
%% Factor AIC
[NumLags, AIC, SCC, HQC,FPE,logL] = VAR_IC(Fr0,12,0);
tableFactorIC = table(AIC,SCC,HQC,FPE,logL);
tableFactorIC.Properties.RowNames = string(1:12);
table2latex(tableFactorIC,'../Tables/FactorIC.tex')
%% Smaller Sample
clearvars -except NumFacs plag PV;
Order = [1,2,3];
FFR1d = 1;
PrepData;
%%
[FAVAR, FAVARopt] = EstimateReducedFAVAR(x2_st, ydata_st, 3, 3, 3, 0, slowcode, 2, namesXY,1,1);
FY_GARCH = FAVAR.fit + FAVAR.residuals;
FY_Gaussian = FAVAR.fit + FAVAR.residualsNonFilter;
%%
t = 25;
J = 400;
[ChowStat, cValue] = IterChow(FY_GARCH,3,3, t, J, DATESnm,1,1,'Chow_Volck_sampleFAVAR_GARCH');
[ChowStat, cValue] = IterChow(FY_Gaussian,3,3, t, J, DATESnm,1,1,'Chow_Volck_sampleFAVAR');
%%
K = 3;               % Number of Factors
M = size(x2, 2);
t1 = size(x2,1);    % time series observations of xdata
p = K+M;             % p is the dimensionality of [Factors, Y]
plag = 3;            % plag is number of lags in the VAR part
X = x2;
Y = ydata_st;
t2 = size(Y,1);    % time series dimension of ydata
%% Cusum-test
name = ["Gaussian";"GARCH"];
for jj = 1:2
    eval(['FY = FY_',char(name(jj)),';']);
    N = size(FY,2);
    NF = 3;
    NL = plag;
    NY = N - NF;
    RHS = FY(NL+1:end,NF+1:N);
    LHS = FY(NL:end-1,:);
    for L = 1:plag-1
        LHS = [LHS, FY(NL-L:end-L-1,:)];
    end
%% 
Testtype = [{'cusum'},{'cusumsq'}];
for j=1:2
        for i=1:3
[~,~,Stat,~,~ ,c1,c2] = cusumtest_modded(LHS, RHS(:,i),'Test',Testtype(j));
Cusum_Stat(:,i,j) = Stat;
Cusum_upper(:,i,j) = c1;
Cusum_lower(:,i,j) = c2;
end
end
%%
index= reshape(1:6, 2, 3).';
Testtype = ["residuals", "squared residuals"];
    figure       
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    namesXY(end,:) = {'Policy Rate'};
    t = 25;
iter = 0;
for j = 1:2
    for i=1:3
        iter = iter+1;
        subplot(3,2,index(iter));   
        plot(DATES(t:end-1,1), Cusum_upper(:,i,j),'LineStyle','--','LineWidth',1.25);
        hold on
        plot(DATES(t:end,1), Cusum_Stat(:,i,j),'LineWidth',1.5);
        hold on 
        plot(DATES(t:end-1,1), Cusum_lower(:,i,j),'LineStyle','--','LineWidth',1.25);
       if j == 1 
           mini = min(min(Cusum_Stat(:,i,j)),min(Cusum_lower(:,i,j)))-20;
           maxi = max(max(Cusum_Stat(:,i,j)),max(Cusum_upper(:,i,j)))+20;
           ylim([mini maxi]);
       else
           mini = min(min(Cusum_Stat(:,i,j)),min(Cusum_lower(:,i,j)))-0.20;
           maxi = max(max(Cusum_Stat(:,i,j)),max(Cusum_upper(:,i,j)))+0.20;
           ylim([mini maxi]);
       end
        title(namesXY(size(slowcode,2)+i),'FontSize',16,'interpreter','latex'); 
       if i == 1
           title({['Cumulative sum of ' Testtype{j}], 'Inflation'},'FontSize',16,'interpreter','latex');
          legend('Critical Value, $5\%$','Test statistic ','interpreter','latex','FontSize',14,'location','Southeast'); 
       else
           title(namesXY(size(slowcode,2)+i),'FontSize',16,'interpreter','latex'); 
       end
       yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
        xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
       betweenShading = [DATES(t:end-1)', fliplr(DATES(t:end-1)')];
       CurvatureShade = [Cusum_lower(:,i,j)', fliplr(Cusum_upper(:,i,j)')];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       recessionplot_RMJ("k",0.20);
       recessionplot_RMJ("k",0.075,"recession",PV);
    end
end
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/Cusum_1983_sample',char(name(jj)),'.pdf']);
end
%%
for jj = 1:2
    eval(['FY = FY_',char(name(jj)),';']);
    N = size(FY,2);
    NF = 3;
    NL = plag;
    NY = N - NF;
    RHS = FY(NL+1:end,1:3);
    LHS = FY(NL:end-1,:);
    for L = 1:plag-1
        LHS = [LHS, FY(NL-L:end-L-1,:)];
    end
%% 
Testtype = [{'cusum'},{'cusumsq'}];
for j=1:2
        for i=1:3
[~,~,Stat,~,~ ,c1,c2] = cusumtest_modded(LHS, RHS(:,i),'Test',Testtype(j));
Cusum_Stat(:,i,j) = Stat;
Cusum_upper(:,i,j) = c1;
Cusum_lower(:,i,j) = c2;
end
end
%%
namesF = ["Factor 1", "Factor 2", "Factor 3"];
index= reshape(1:6, 2, 3).';
Testtype = ["residuals", "squared residuals"];
    figure       
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    namesXY(end,:) = {'Policy Rate'};
    t = 25;
iter = 0;
for j = 1:2
    for i=1:3
        iter = iter+1;
        subplot(3,2,index(iter));   
        plot(DATES(t:end-1,1), Cusum_upper(:,i,j),'LineStyle','--','LineWidth',1.25);
        hold on
        plot(DATES(t:end,1), Cusum_Stat(:,i,j),'LineWidth',1.5);
        hold on 
        plot(DATES(t:end-1,1), Cusum_lower(:,i,j),'LineStyle','--','LineWidth',1.25);
       if j == 1 
           mini = min(min(Cusum_Stat(:,i,j)),min(Cusum_lower(:,i,j)))-20;
           maxi = max(max(Cusum_Stat(:,i,j)),max(Cusum_upper(:,i,j)))+20;
           ylim([mini maxi]);
       else
           mini = min(min(Cusum_Stat(:,i,j)),min(Cusum_lower(:,i,j)))-0.20;
           maxi = max(max(Cusum_Stat(:,i,j)),max(Cusum_upper(:,i,j)))+0.20;
           ylim([mini maxi]);
       end
        title(namesF(i),'FontSize',16,'interpreter','latex'); 
       if i == 1
           title({['Cumulative sum of ' Testtype{j}], 'Factor 1'},'FontSize',16,'interpreter','latex');
          legend('Critical Value, $5\%$','Test statistic ','interpreter','latex','FontSize',14,'location','Southeast'); 
       else
           title(namesF(i),'FontSize',16,'interpreter','latex'); 
       end
       yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
        xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
       betweenShading = [DATES(t:end-1)', fliplr(DATES(t:end-1)')];
       CurvatureShade = [Cusum_lower(:,i,j)', fliplr(Cusum_upper(:,i,j)')];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       recessionplot_RMJ("k",0.20);
       recessionplot_RMJ("k",0.075,"recession",PV);
    end
end
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/CusumFactors_1983_sample',char(name(jj)),'.pdf']);
end