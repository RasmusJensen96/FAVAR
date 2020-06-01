function PlotResidualFAVARGarch(FAVAR,FAVARopt,DATES,Name)

epsilon = FAVAR.residualsNonFilter;
sigma   = FAVAR.volatility;
z       = FAVAR.residuals;

L = FAVAR.nlag;
N = FAVAR.nvar;
Names = FAVARopt.vnames;
DATES = DATES(L+1:end);

PV = [datenum(1979,08,06), datenum(1987,08,11)];
figure;
foo = {[1; 3]; 5; [7; 9]; 11; [13; 15]; 17; [2; 4]; 6; [8; 10]; 12; [14; 16]; 18};
set(0,'DefaultAxesColorOrder','remove');
indicator = 1;
for i=1:N
            subplot(9,2,foo{indicator})
                x = plot(DATES, epsilon(:,i));
                hold on
                y = plot(DATES, sigma(:,i),'LineWidth',1.25);
                title(Names(i), 'interpreter', 'latex', 'FontSize',14);
                   set(gca,'xtick',[])
                   set(gca,'xticklabel',[])
                   ylabel('Residuals','Interpreter','latex','FontSize',14);
                   if i == 1
                   legend([y,x],'$\sigma_t$','$\epsilon_t=\sigma_t z_t$','Location','best','interpreter','latex');
                   end
                   ylim([min(min(epsilon(:,i)), min(epsilon(:,i)))-0.5*std(epsilon(:,i)), max(max(epsilon(:,i)), max(sigma(:,i)))+0.5*std(sigma(:,i))]);
                   recessionplot_RMJ("k",0.20);
                   recessionplot_RMJ("k",0.075,"recession",PV); 
                   yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
                   xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
                indicator = indicator + 1 ;
            subplot(9,2,foo{indicator})
                plot(DATES, z(:,i),'k');
                recessionplot_RMJ("k",0.20);
                recessionplot_RMJ("k",0.075,"recession",PV); 
               % ylabel('Residuals');
               if indicator == 6 | indicator == 12;
                    xlabel = ('Date');
               end
               yline(0,'--', 'LineWidth',1.25);
               yAX = get(gca,'YAxis');yAX.FontSize = 12;yAX.TickLabelInterpreter = "latex";
               xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
               ylabel({'$z_t$'},'Interpreter','latex');
                if i ~= [3, 6]
                   set(gca,'xtick',[])
                   set(gca,'xticklabel',[])
               end
                indicator = indicator +1;
end
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/', Name,'.pdf']);