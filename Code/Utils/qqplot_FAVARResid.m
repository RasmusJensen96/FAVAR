function  qqplot_FAVARResid(FAVAR, FAVARopt,name)
N = size(FAVARopt.vnames, 1);
figure;
for i = 1:N
    subplot(2, 3, i)
    Y = qqplot_Mod(FAVAR.residuals(:,i));
    title(FAVARopt.vnames{i},'interpreter','latex','FontSize',30);
    yAX = get(gca,'YAxis');yAX.FontSize = 30;yAX.TickLabelInterpreter = "latex";
    xAX = get(gca,'XAxis');xAX.FontSize = 30;xAX.TickLabelInterpreter = "latex";
    if i == 1 | i==4
        ylabel('Empirical Quantiles');
    else
        ylabel(' ');
    end
    %if i == 4 | i==5 | i == 6;
    %    xlabel('Standard Normal Quantiles');
    %else
        xlabel(' ');
    %end
end
% get the original size of figure before the legends are added
%figuresize =  get(gcf, 'position');

% set unit for figure size to inches
set(gcf, 'unit', 'inches');
% get the original size of figure before the legends are added
%figure_size =  get(gcf, 'position');
a = annotation('textbox',[0.1 0.1 1.0 0.23],'string',{"Standard Normal Quantiles"},'Interpreter','latex','FontSize',30,'EdgeColor','none');
a.Units = 'inches'; a.Position = [3 0.1 5.25 .3];


    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
    saveas(gcf,['../Figures/',name,'.pdf']);
