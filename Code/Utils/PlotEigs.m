function [] = PlotEigs(FAVAR,name)
figure % Plot eigenvals
plot(FAVAR.Eigs,'o');
ylabel('Imaginary','Interpreter','latex','fontsize',14); xlabel('Real','Interpreter','latex','fontsize',14)
xline(0,'--');yline(0,'--');
hold on
 pos = [-1 -1 2 2];
   rectangle('Position',pos,'Curvature',[1 1])
   rectangle('Position',[-.75 -.75 1.5 1.5],'Curvature',[1 1])
   rectangle('Position',[-.5 -.5 1 1],'Curvature',[1 1])
   rectangle('Position',[-.25 -.25 .5 .5],'Curvature',[1 1])
   axis equal
   title('Eigenvalues of the FAVAR','Interpreter','latex','FontSize',16);
   yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
   xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
set(gcf, 'PaperPosition', [0 0 40 30]);
set(gcf, 'PaperSize', [40 30]);
saveas(gcf,['../Figures/',name,'.pdf']);
end