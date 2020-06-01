function FAVAR_Resp(FAVAR, FAVARopt, IRFMED_Chol, IRFINF_Chol, IRFSUP_Chol,FileName, FFR1d, factors)

if ~exist('FAVARopt.standardized','var');
    FAVARopt.standardized = 1;
end
if FAVARopt.standardized == true
   foo_str = {'Deviation from mean, std. units'};
else
   foo_str = {'Deviation from mean'}; 
end

vnames = FAVARopt.vnames;
K = FAVAR.NumberFactors;

% entFFR     = find([vnames{:}] == "FFR");
% entIP      = find([vnames{:}] == "IP");
% entInf     = find([vnames{:}] == "Inflation");

entFFR     = find(strcmp(vnames, 'FFR'));
entIP      = find(strcmp(vnames, 'IP'));
entInf     = find(strcmp(vnames, 'CPI'));
entfoo     = [entIP, entInf];

if string(FAVARopt.ident) == "oir"
    entplotindex = [entInf, entIP, entFFR];
else
    entplotindex = [entIP, entFFR, entInf];
end
if FFR1d == 1
    IRFMED_Chol(:,entFFR,:) = cumsum(IRFMED_Chol(:,entFFR,:));
    IRFINF_Chol(:,entFFR,:) = cumsum(IRFINF_Chol(:,entFFR,:));
    IRFSUP_Chol(:,entFFR,:) = cumsum(IRFSUP_Chol(:,entFFR,:));
end

IRFMED_Chol(:,entfoo,:) = exp(cumsum(IRFMED_Chol(:,entfoo,:)))-1;
IRFINF_Chol(:,entfoo,:) = exp(cumsum(IRFINF_Chol(:,entfoo,:)))-1;
IRFSUP_Chol(:,entfoo,:) = exp(cumsum(IRFSUP_Chol(:,entfoo,:)))-1;

% IRFMED_Chol(:,entfoo,:) = cumsum(IRFMED_Chol(:,entfoo,:));
% IRFINF_Chol(:,entfoo,:) = cumsum(IRFINF_Chol(:,entfoo,:));
% IRFSUP_Chol(:,entfoo,:) = cumsum(IRFSUP_Chol(:,entfoo,:));

%IRFMED_Chol(:,K+1:end-1,:) = exp(cumsum(IRFMED_Chol(:,K+1:end-1,:)))-1;
%IRFINF_Chol(:,K+1:end-1,:) = exp(cumsum(IRFINF_Chol(:,K+1:end-1,:)))-1;
%IRFSUP_Chol(:,K+1:end-1,:) = exp(cumsum(IRFSUP_Chol(:,K+1:end-1,:)))-1;

if FAVARopt.CIMethod
    textstr = ['Confidence interval: $\pm',num2str(FAVARopt.nStd),'\sigma$'];
else
    textstr = ['Confidence interval: $',num2str(FAVARopt.pctg),'\%$'];
end
if isempty(FAVARopt.snames)
    snames = FAVARopt.vnames;
else
    snames = FAVARopt.snames;
end
ynames = FAVARopt.vnames(entplotindex);
%% IRF plot chol;
if factors == 1
NumShock = size(FAVAR.ENDO,2);
numResp = (NumShock-K);
totalSup = numResp * NumShock;
foo = reshape(1:totalSup, NumShock, numResp).';
figure       
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    iter = 0;
for k = 1:NumShock
     for i=1:numResp
        iter = iter + 1; 
        subplot(size(ynames, 1),size(snames,1),foo(iter));   
         yt =  plot(IRFMED_Chol(:,entplotindex(i),k),'LineWidth',1.5);
               title(['Shock: ', snames{k}],'Interpreter','latex','FontSize',16);
         hold on
              plot(IRFINF_Chol(:,entplotindex(i),k),'--','LineWidth',1.25);
         hold on
            xt =  plot(IRFSUP_Chol(:,entplotindex(i),k),'--','LineWidth',1.25);
              yline(0,'-','LineWidth',1);
              yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
              xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
       if i == 3;
          xlabel('Periods','interpreter','latex','FontSize',14);
       end
       if foo(iter) == 1 | foo(iter) == 7 | foo(iter) == 13
           ylabel([{['Response: ',  ynames{i}]}; foo_str],'interpreter','latex','FontSize',14);
       end
       betweenShading = [1:FAVARopt.nsteps fliplr(1:FAVARopt.nsteps)];
       CurvatureShade = [IRFINF_Chol(:,entplotindex(i),k)', fliplr(IRFSUP_Chol(:,entplotindex(i),k)')];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       xlim([1 FAVARopt.nsteps]);
end
end
% set unit for figure size to inches
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
    
   saveas(gcf,['../Figures/IRFs',FileName,'wF.pdf']);
else
    
NumShock = size(FAVAR.ENDO,2)-K;
numResp = NumShock;
totalSup = numResp * NumShock;
    ynames = FAVARopt.vnames(entplotindex);
    if size(ynames,1) < size(ynames,2)
           ynames = ynames'; 
           snames = snames';
    end
    snames = snames(entplotindex);
    foo = reshape(1:totalSup, NumShock, numResp).';
figure       
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    iter = 0;
for k = 1:NumShock
     for i=1:numResp
        iter = iter + 1; 
        subplot(size(ynames, 1),size(snames,1),foo(iter));   
         yt =  plot(IRFMED_Chol(:,entplotindex(i),entplotindex(k)),'LineWidth',1.5);
               title(['Shock: ', snames{k}],'Interpreter','latex','FontSize',16);
         hold on
              plot(IRFINF_Chol(:,entplotindex(i),entplotindex(k)),'--','LineWidth',1.25);
         hold on
            xt =  plot(IRFSUP_Chol(:,entplotindex(i),entplotindex(k)),'--','LineWidth',1.25);
              yline(0,'-','LineWidth',1);
              yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
              xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
       if i == 3;
          xlabel('Periods','interpreter','latex','FontSize',14);
       end
       if foo(iter) == 1 | foo(iter) == 1+NumShock | foo(iter) == 1 + 2 * NumShock
           ylabel([{['Response: ',  ynames{i}]}; foo_str],'interpreter','latex','FontSize',14);
       end
       betweenShading = [1:FAVARopt.nsteps fliplr(1:FAVARopt.nsteps)];
       CurvatureShade = [IRFINF_Chol(:,entplotindex(i),entplotindex(k))', fliplr(IRFSUP_Chol(:,entplotindex(i),entplotindex(k))')];
       fill(betweenShading, CurvatureShade,'k','FaceAlpha',0.05,'HandleVisibility','off');
       xlim([1 FAVARopt.nsteps]);
end
end
% set unit for figure size to inches
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
    
 saveas(gcf,['../Figures/IRFs',FileName,'.pdf']);
end