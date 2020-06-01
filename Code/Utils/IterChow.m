function [ChowStat, cValue] = IterChow(FY, NumFacs,NumLags,StartObs,EndObs,yearlab,SavePlot,RecBands,name)
% Conducts an iterative Chow-test for the VAR-part of the FASVAR:
% Input = FY: stacked vector of time-seris in the FASVAR:
%         NumFacs: Number of factors
%         NumLags: Number of lags

N = size(FY,2);
NF = NumFacs;
NL = NumLags;
NY = N - NF;

RHS = FY(NL+1:end,NF+1:N);
LHS = FY(NL:end-1,:);

for L = 1:NumLags-1
    LHS = [LHS, FY(NL-L:end-L-1,:)];
end

for series=1:NY
    iter = 0;
    for bp=StartObs:EndObs
        iter = iter + 1;
        [~,~,stat,cValue] = chowtest(LHS, RHS(:,series), bp,'alpha',0.05,'Intercept',false);
        ChowStat(iter,series) = stat;
    end
end

t = datetime(yearlab(StartObs:EndObs),'ConvertFrom','datenum');
t = t + 365.25 * 1900;
names = [{'Inflation'},{'Industrial Production'},{'Policy Rate'}];

if RecBands == true
% Auxiliary to add Volcker regime:
Inaug = datenum(1979,08,06);
%Resgn = datetime(1987,08,11);
Resgn = datenum(1987,08,11);
PV = [Inaug Resgn];

end
figure;
for i=1:NY
subplot(NY,1,i)
plot(t, ChowStat(:,i),'LineWidth',1.25);
ylim([min(ChowStat(:,i))-std(ChowStat(:,i)), max(ChowStat(:,i))+std(ChowStat(:,i))]);
title(names{i},'interpreter','latex','FontSize',16);
yline(cValue, '--', 'LineWidth',1.25);
ylabel("Chow-Statistic",'FontSize',14,'interpreter','latex');
if RecBands == true
    recessionplot_RMJ("k",0.20);
        %recessionplot_RMJ([0 0.4470 0.7410],0.2,"recession",PV);
    recessionplot_RMJ("k",0.075,"recession",PV); 
end
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
if i == 1
    legend('Chow-Statistic', 'Critical Value, $\alpha = 0.05$','interpreter','latex','FontSize',12,'Location','northeast');
end
end
xlabel("Date",'FontSize',14,'interpreter','latex');
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);

if SavePlot == true
saveas(gcf,['../Figures/', name, '.pdf'])
end

end