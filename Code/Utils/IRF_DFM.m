function IRF_DFM(FAVAR, FAVARopt, MEDx, INFx, SUPx,  tcodes, var_numbers, var_names, prefix, Order, FFR1d)
% Function plotting the responses of the factor panel to FAVAR innovations
% Rasmus M. Jensen 2020

var_names = strrep(var_names,'&','\&');
var_names = strrep(var_names,'\\&','\&');

if FAVARopt.CIMethod
    textstr = ['Confidence interval: $\pm',num2str(FAVARopt.nStd),'\sigma$'];
else
    textstr = ['Confidence interval: $',num2str(FAVARopt.pctg),'\%$'];
end


if isempty(FAVARopt.snames)
    vnames = FAVARopt.vnames;
else
    vnames = FAVARopt.snames;
end

for i = 1:3
    if Order(i) == 1;
        tcy(i) = 5;
    elseif Order(i) == 2;
        tcy(i) = 5;
    else
        if FFR1d == 1;
        tcy(i) = 2;
        else
            tcy(i) = 1;
        end
    end
end


nlag = num2str(FAVAR.nlag);
tcodes = [tcodes, tcy];

for j=1:size(MEDx,3)
    figure
name = vnames{j};
%% DFM Cholesky
DFM = MEDx(:,:,j);
DFM_Inf = INFx(:,:,j);
DFM_Sup = SUPx(:,:,j);
[IRFs, Upper, lower] = DetransformIRF(DFM', DFM_Sup', DFM_Inf', tcodes);
ImpRepXD = IRFs';
ConfInfD = lower';
ConfSupD = Upper';

for k = 1:2
    set(0,'DefaultAxesColorOrder',[0 0 0]);
% IRF factors:
    for i=1:9
        subplot(3,3,i)   
        plot(ConfInfD(var_numbers(i),:),'LineStyle','--','LineWidth',1.25);
        hold on
        yt = plot(ImpRepXD(var_numbers(i),:),'LineWidth',1.5);
        hold on 
        xt = plot(ConfSupD(var_numbers(i),:),'LineStyle','--','LineWidth',1.25);
        hold on
        yline(0,'-','LineWidth',1);
             yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
             xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
        title([var_names(i) ' to shock to' name],'FontSize',16,'interpreter','latex');  
       if i == 7 | i == 8 | i == 9
          xlabel('Periods','interpreter','latex','FontSize',14);
       end
       if i == 1 | i == 4 | i == 7
           ylabel('Deviation from mean, std. units','interpreter','latex','FontSize',14);
       end
      % if i == 3
       %   legend('Median Response','CI-Bounds $68\%$','interpreter','latex','FontSize',14,'location','best'); 
       % end
       betweenShading = [1:size(ConfInfD(var_numbers(i),:),2) fliplr(1:size(ConfInfD(var_numbers(i),:),2))];
       CurvatureShade = [ConfInfD(var_numbers(i),:), fliplr(ConfSupD(var_numbers(i),:))];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       xlim([1, size(ConfInfD(var_numbers(i),:),2)]);
    end
    
set(gcf, 'unit', 'inches');
% get the original size of figure before the legends are added
%figure_size =  get(gcf, 'position');
lg = legend([yt, xt], 'Median Response',textstr,'interpreter','latex','FontSize',14,'location','southoutside','Orientation','horizontal'); 

set(lg,'Units','inches',...
    'Position',[3.5 0.1 1.0 0.23],...
    'Interpreter','latex',...
    'FontSize',14);
    
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
saveas(gcf,strrep(join(string(['../Figures/IRFs',prefix,'DFMshockto',name,num2str(nlag),'Lags.pdf'])),' ',''));
end
end
end