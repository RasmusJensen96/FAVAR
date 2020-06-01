function VARhdplot_GP(HD,VARopt,dates)
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

% modified by Giovanni Pellegrino

%% Check inputs
%===============================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
% If there is VARopt check that vnames is not empty
vnames = VARopt.vnames;
if isempty(vnames)
    error('You need to add label for endogenous variables in VARopt');
end


%% Check inputs Define some parameters
%===============================================
filename = [VARopt.figname 'HD_'];
quality = VARopt.quality;
suptitle = VARopt.suptitle;
pick = VARopt.pickHD;

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

HD_pickShock=VARopt.HD_pickShock  ;

%% Plot
%===============================================
% FigSize
for ii=pick:nvars                
    figure
    colormap(parula)
    bar(dates,HD.shock(:,HD_pickShock,ii));hold on
    plot(dates,HD.endo(:,pick)-HD.init(:,pick)-HD.const(:,pick)-HD.trend(:,pick)-HD.trend2(:,pick)-HD.exo(:,pick))
%     xlim([1 nsteps]);
    title([vnames{ii} ' (stoch.comp.)'], 'FontWeight','bold','FontSize',10); 
    % Save
    FigName = [filename num2str(ii)];
    if quality 
        if suptitle==1
            Alphabet = char('a'+(1:nvars)-1);
            SupTitle([Alphabet(jj) ') HD of '  vnames{ii}])
        end
        opt = LegOption; LegSubplot(vnames,opt);
        set(gcf, 'Color', 'w');
        export_fig(FigName,'-pdf','-png','-painters')
    else
        if suptitle==1
            Alphabet = char('a'+(1:nvars)-1);
            SupTitle([' HD of '  vnames{ii}])
        end
        legend([vnames(HD_pickShock)])
        print('-dpng','-r100',FigName);
        print('-dpdf','-r100',FigName);
    end
%     clf('reset');
end

% close all