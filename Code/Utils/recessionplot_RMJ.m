function varargout = recessionplot(colourchosen,alphachosen,varargin)
% Modified 04-02-2020 by R. M. Jensen
%RECESSIONPLOT Add recession bands to time series plot
%
% Syntax:
%
%   recessionplot
%   recessionplot(param,val,...)
%   hBands = recessionplot(...)
%
% Description:
%
%   Overlay shaded recession bands on a time series plot.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME            VALUE
%
%   'axes'          Handle to axes displaying a time series plot with
%                   serial date numbers on the horizontal axis. The default
%                   is gca.
%
%	'recessions'	numRecessions-by-2 matrix of serial date numbers
%                   indicating the beginning (first column) and end (second
%                   column) of historical recessions. The default is the
%                   U.S. recession data in Data_Recessions.mat, reported
%                   by the National Bureau of Economic Research.
%
% Output Arguments:
%
%   hBands - Vector of handles to the recession bands.
%
% Notes:
%
%   o RECESSIONPLOT requires dates on the horizontal axis of a time series
%     plot to be expressed as serial date numbers. To convert other date
%     information to this format before plotting, use the DATENUM command.
%
%   o Use the output handles to change the color and opacity of the
%     recession bands by setting their 'FaceColor' and 'FaceAlpha'
%     properties. This may be necessary to achieve satisfactory displays
%     when working with certain monitors and projectors.
%
% Example:
%
%   load Data_CreditDefaults   
%   X0 = Data(:,1:4);
%   T0 = size(X0,1);
%   
%   % Convert dates to serial date numbers:
% 
%   dates = datenum([dates,ones(T0,2)]);
% 
%   % Create plot:
%   
%   fig = figure;
%   ax = axes(fig);
%   plot(ax,dates,X0,'LineWidth',2)
%   set(ax,'XTick',dates(1:2:end))
%   datetick(ax,'x','yyyy','keepticks')
%   xlabel(ax,'Year') 
%   ylabel(ax,'Level')
%   axis(ax,'tight')
%  
%   % Add recession bands:
% 
%   recessionplot('axes',ax)
%
% References:
% 
%   [1] National Bureau of Economic Research. "US Business Cycle Expansions
%       and Contractions." http://www.nber.org/cycles.html.
%
% See also DATENUM.

% Copyright 2012-2017 The MathWorks, Inc.

% Parse inputs and set defaults:

% Modified by Rasmus M. Jensen 2020

nargoutchk(0,1);

parseObj = inputParser;
parseObj.addParameter('axes',get(get(0,'CurrentFigure'),'CurrentAxes'),@axesCheck);
parseObj.addParameter('recessions',[],@recessionsCheck);

parseObj.parse(varargin{:});

hAx = parseObj.Results.axes;
Recessions = parseObj.Results.recessions;

if isempty(hAx)
        
	error(message('econ:recessionplot:NoCurrentAxes'))

end

if isempty(colourchosen)
    colourchosen = "k";
end

if isempty(alphachosen)
    alphachosen = 0.1;
end

if isempty(Recessions)
    
    load Data_Recessions
    
end

dateRange = get(hAx,'XLim'); % Date range of current axes
dataRange = get(hAx,'YLim'); % Data range of current axes

if isdatetime(dateRange)
    if isnumeric(Recessions)
        Recessions = datetime(Recessions, 'ConvertFrom', 'datenum');
    end
else
    if isdatetime(Recessions)
        Recessions = datenum(Recessions);
    end
end

inPlot = (Recessions(:,1) < dateRange(2)) & (Recessions(:,2) > dateRange(1));
RPlot = Recessions(inPlot,:); % Recessions to plot

numRecessions = size(RPlot,1);
hBands = gobjects(numRecessions,1);  % Recession band handles
nextPlot = get(hAx, 'NextPlot');     % Save the current NextPlot property
set(hAx, 'NextPlot', 'add');         % Ensure recession patches are added

for n = 1:numRecessions

       hBands(n) = fill(hAx,[RPlot(n,1),RPlot(n,1),RPlot(n,2),RPlot(n,2)], ...
                            [dataRange,fliplr(dataRange)],colourchosen,'FaceAlpha',alphachosen,'EdgeColor','none','HandleVisibility','off');

end

set(hAx, 'NextPlot', nextPlot);      % Restore the original NextPlot property
xlim(hAx,dateRange)

if nargout > 0
    
    varargout{1} = hBands;
    
end

%-------------------------------------------------------------------------
% Check value of 'axes' parameter
function OK = axesCheck(hAx)

if ~ishghandle(hAx,'axes')
    
    error(message('econ:recessionplot:AxesHandleInvalid'))

else
    
    OK = true;
    
end

%-------------------------------------------------------------------------
% Check value of 'recessions' parameter
function OK = recessionsCheck(Recessions)
    
if (~isnumeric(Recessions)&&~isdatetime(Recessions)) || ~ismatrix(Recessions) || size(Recessions,2) ~= 2

    error(message('econ:recessionplot:RecessionsMatrixInvalid'))

else

    OK = true;

end