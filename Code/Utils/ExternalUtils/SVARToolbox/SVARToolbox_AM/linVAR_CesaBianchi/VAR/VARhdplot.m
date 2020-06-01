function VARhdplot(HD,VARopt, FAVAR, DATES, Modelname)
% =======================================================================
% Plot the HD shocks computed with VARhd
% =======================================================================
% VARhdplot(HD,VARopt)
% -----------------------------------------------------------------------
% INPUT
%   - HD: structure from VARhd
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

Endo = FAVAR.ENDO;
%% Check inputs
%===============================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
% If there is VARopt check that vnames is not empty
vnames = string(VARopt.vnames);
if isempty(vnames)
    error('You need to add label for endogenous variables in VARopt');
end


%% Check inputs Define some parameters
%===============================================
pick = VARopt.pick;

% Initialize HD matrix
[nsteps, nvars, nshocks] = size(HD.shock);

% If one shock is chosen, set the right value for nshocks
if pick<0 || pick>nvars
    error('The selected shock is non valid')
else
    if pick==0
        pick=1;
    else
        nshocks = pick;
    end
end

nshocks = 6
%% Plot
%===============================================
% FigSize
for ii=pick:nvars        
    name = join(['HD_',vnames(ii),'_Model',Modelname],'');
    figure
    %colormap(parula)
    H = BarPlot(HD.shock(:,1:nshocks,ii), DATES);
    hold on;
    Y       = plot(DATES, Endo(:,ii),'k','LineWidth',1,'DisplayName','Observation');
    HJ      = H(1,:)';
    HJ(7,:) = Y;
    title([vnames(ii)], 'Interpreter','latex','FontSize',14); 
    yAX = get(gca,'YAxis');yAX.FontSize = 14;yAX.TickLabelInterpreter = "latex";
    xAX = get(gca,'XAxis');xAX.FontSize = 14;xAX.TickLabelInterpreter = "latex";
    legend(HJ,[vnames; 'Observation'], 'interpreter','latex','Location','bestoutside','Orientation','horizontal');
    % Save
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
    recessionplot_RMJ("k",0.20);
    saveas(gcf,join(['../Figures/HistoricalDecomp/', name ,'.pdf'],''));
%     clf('reset');
end
%hold on
%plot(DATES, Endo(:,pick),'k','LineWidth',1.5);
%legend([vnames; {'Observed'}], 'interpreter','latex','Location','best');

% close all
