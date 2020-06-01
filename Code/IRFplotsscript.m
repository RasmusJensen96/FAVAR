set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
set(0,'DefaultAxesColorOrder',[0 0 0]);

Pre = load('IRFMedbqPre');
Post = load('IRFMedbqPost');
Full = load('IRFMedbqFull');
Bernanke = load('IRFMedbqBern');

figure;
plot(Post.IRFMED_bq(:,6,3),'LineWidth',1.5)
hold on
plot(Pre.IRFMEDPre_bq(:,6,3),'--','LineWidth',1.5)
plot(Full.IRFMED_bq(:,6,3),'.-','LineWidth',1.5,'MarkerSize',15)
plot(Bernanke.IRFMED_bq(1:24,6,3),':','LineWidth',2)
yline(0,'-')
legend(["1984:2019","1959:1984", "1959:2019", "1959:2001"],'interpreter','latex','FontSize',24); 
yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/IRFsSpendingRatediff.pdf']);
figure
plot(exp(cumsum(Pre.IRFMEDPre_bq(:,6,3)))-1,'--','LineWidth',1.5)
hold on 
plot(exp(cumsum(Post.IRFMED_bq(:,6,3)))-1,'LineWidth',1.5)
plot(exp(cumsum(Full.IRFMED_bq(:,6,3)))-1,'.-','LineWidth',1.5,'MarkerSize',15)
plot(exp(cumsum(Bernanke.IRFMED_bq(1:24,6,3)))-1,':','LineWidth',2)
yline(0,'-')
yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/IRFsSpendingRatediffCumum.pdf']);
% %%
% t = 24;
% set(0,'DefaultAxesColorOrder',[0 0 0]);
% Pre = load('IRFMedCholPre');
% Post = load('IRFMedCholPost');
% Full = load('IRFMedCholFull');
% Bernanke = load('IRFMedCholBern');
% %%
% figure;
% plot(Post.IRFMED_Chol(1:t,4,6),'LineWidth',1.5)
% hold on
% plot(Pre.IRFMED_Chol(1:t,4,6),'--','LineWidth',1.5)
% plot(Full.IRFMED_Chol(1:t,4,6),'.-','LineWidth',1.5,'MarkerSize',15)
% plot(Bernanke.IRFMED_Chol(1:t,4,6),':','LineWidth',2)
% yline(0,'-')
% legend(["1984:2019","1959:1984", "1959:2019", "1959:2001"],'interpreter','latex','FontSize',24); 
% yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
% xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
% set(gcf, 'PaperPosition', [0 0 40 30]);
% set(gcf, 'PaperSize', [40 30]);
% saveas(gcf,['../Figures/IRFsCholInfl2FFRdiff.pdf']);
% 
% figure
% plot(exp(cumsum(Post.IRFMED_Chol(1:t,4,6)))-1,'LineWidth',1.5)
% hold on
% plot(exp(cumsum(Pre.IRFMED_Chol(1:t,4,6)))-1,'--','LineWidth',1.5)
% plot(exp(cumsum(Full.IRFMED_Chol(1:t,4,6)))-1,'.-','LineWidth',1.5,'MarkerSize',15)
% plot(exp(cumsum(Bernanke.IRFMED_Chol(1:t,4,6)))-1,':','LineWidth',2)
% yline(0,'-')
% yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
% xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
% set(gcf, 'PaperPosition', [0 0 40 30]);
% set(gcf, 'PaperSize', [40 30]);
% saveas(gcf,['../Figures/IRFsCholCPI2FFRdiff.pdf']);