function [] = FAVAR_Cusum(FAVAR,DATES,namesXY,K,figname,Chol)

plag = FAVAR.nlag;
PV = [datenum(1979,08,06), datenum(1987,08,11)]; % Paul Volcker office period.

FY = FAVAR.fit + FAVAR.residuals;

N = size(FY,2);
%RHS = FY(plag+1:end,K+1:N);
if Chol == 0
    RHS = FY(plag+1:end,K+1:end);
else
    RHS = FY(plag+1:end,[2,3,6]);
end
%RHS = FY(plag+1:end,1:3);
LHS = FY(plag:end-1,:);
for L = 1:plag-1
    LHS = [LHS, FY(plag-L:end-L-1,:)];
end

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
Sz = size(DATES,1) - size(Cusum_Stat,1);

set(0,'DefaultAxesColorOrder',[0 0 0]);
index= reshape(1:6, 2, 3)';
Testtype = ["residuals", "squared residuals"];
    figure       
    t = 22;
iter = 0;
for j = 1:2
    for i=1:3
        iter = iter+1;
        subplot(3,2,index(iter));   
        plot(DATES(Sz+1:end-1,1), Cusum_upper(:,i,j),'LineStyle','--','LineWidth',1.25);
        hold on
        plot(DATES(Sz+1:end,1), Cusum_Stat(:,i,j),'LineWidth',1.5);
        hold on 
        plot(DATES(Sz+1:end-1,1), Cusum_lower(:,i,j),'LineStyle','--','LineWidth',1.25);
        %title(namesXY(size(namesX,1)+i),'FontSize',16,'interpreter','latex'); 
        title(namesXY(size(FAVAR.DFM.SIGMA,1)-3+i),'FontSize',16,'interpreter','latex'); 
       if i == 1
          title({['Cumulative sum of ' Testtype{j}], 'Inflation'},'FontSize',16,'interpreter','latex');
          legend('Critical Value, $5\%$','Test statistic ','interpreter','latex','FontSize',14,'location','southeast'); 
       else
           title(namesXY(size(FAVAR.DFM.SIGMA,1)-3+i),'FontSize',16,'interpreter','latex'); 
       end
       betweenShading = [DATES(Sz+1:end-1)', fliplr(DATES(Sz+1:end-1)')];
       CurvatureShade = [Cusum_lower(:,i,j)', fliplr(Cusum_upper(:,i,j)')];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       %recessionplot;
       recessionplot_RMJ("k",0.20);
       recessionplot_RMJ("k",0.075,"recession",PV); 
       yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
       xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
       if j == 2
       ylim([-0.1268 1]); % Recessionband correction when auto-exporting
       end
    end
end
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/',figname,'.pdf']);