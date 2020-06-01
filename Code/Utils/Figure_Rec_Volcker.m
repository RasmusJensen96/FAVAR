clear all
%%%
Chairman  = ["Paul Volcker"];
%Inaug = datetime(1979,08,06);
Inaug = datenum(1979,08,06);
%Resgn = datetime(1987,08,11);
Resgn = datenum(1987,08,11);
PV = [Inaug Resgn];
%%
ydata = readmatrix("VAR_Data.xlsx");
ydata = ydata(1:end-3,1:3);
yearlab = readmatrix("Dates.xlsx");
t = datetime(yearlab,'ConvertFrom','datenum');
t = t + 365.25 * 1900;
%% 
load Data_Recessions;
%%
VolckerRec2 = Recessions(end-3,:);
%%
set(0,'DefaultAxesColorOrder','factory');
y = plot(t, (ydata(:,1:3)-mean(ydata(:,1:3)))./std(ydata(:,1:3)),'LineWidth',1.25);
%y = plot(t, ydata,'LineWidth',1.5);
recessionplot_RMJ("k",0.20);
%recessionplot_RMJ([0 0.4470 0.7410],0.2,"recession",PV);
recessionplot_RMJ("k",0.075,"recession",PV);
y1 = xline(datetime(max(VolckerRec2),'ConvertFrom','datenum'),'--','LineWidth',1.25);
ylabel("Standard deviation units",'FontSize',14,'interpreter','latex');
xlabel("Date",'FontSize',14,'interpreter','latex');
legend([y; y1],["Inflation", "Industrial Production", "Policy-Rate", "Structural Break"]','interpreter','latex','FontSize',12)
yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,'../Figures/Structural_Break.pdf')