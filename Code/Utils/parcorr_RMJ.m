function varargout = parcorr_RMJ(varargin)
%PARCORR Sample partial autocorrelation
%
% Syntax:
%
%   [pacf,lags,bounds] = parcorr(y)
%   [pacf,lags,bounds] = parcorr(y,param,val,...)
%   [pacf,lags,bounds,h] = parcorr(...)
%   [pacf,lags,bounds,h] = parcorr(ax,...)
%   parcorr(...)
%
% Description:
%
%   Compute the sample partial autocorrelation function (PACF) of a
%   univariate, stochastic time series y. PARCORR optionally plots the PACF
%   sequence with confidence bounds.
%
% Input Arguments:
%
%   y - Vector of observations of a univariate time series. The last
%       element of y contains the most recent observation.
%
%   ax - Axes object in which to plot. If unspecified, PARCORR plots to the
%       current axes (gca).
%
% Optional Input Parameter Name/Value Pairs:
%
%  'NumLags' Positive integer that determines the number of lags at which the
%            PACF is computed. The lags used to compute the PACF are 0:NumLags.
%            The default is min[20,N-1], where N is the effective sample size
%            of y
%
%  'NumAR'   For computing confidence bounds, a nonnegative integer less
%            than NumLags specifying the number of lags in a theoretical
%            AR(NumAR) model of y. For lags > NumAR, PARCORR computes the
%            standard error using 1/sqrt(N), for Gaussian white noise,
%            under the model assumption. The default is 0.
%
%  'NumSTD'  For computing confidence bounds, a nonnegative scalar multiple
%            specifying an interval of +/-(NumSTD) times the computed
%            standard error. The default is 2 (approximate 95% confidence).
%
%  'Method'  String or character vector indicating the method by which the
%            PACF is computed. Values are 'OLS', to use ordinary least
%            squares, and 'Yule-Walker', to use the Yule-Walker equations.
%            If the series y is fully observed, the default is 'OLS';
%            otherwise the default is 'Yule-Walker'.
%
% Output Arguments:
%
%   pacf - Sample PACF. Vector of length NumLags+1 of values computed at lags
%       0,1,2,...,NumLags. For all y, pacf(1) = 1 at lag 0.
%
%   lags - Vector of lag numbers of length NumLags+1 used to compute pacf.
%
%   bounds - Two-element vector of approximate upper and lower confidence
%       bounds, assuming that y is an AR(NumAR) process.
%
%   h - Vector of handles to plotted graphics objects. PARCORR plots the
%       PACF when the number of output arguments is 0 or 4.
%
% Notes:
%
%  o Specify missing observations of y using NaN. PARCORR treats these
%    values as "missing completely at random."
%
%  o OLS estimates of PACF coefficients are available only if y is fully
%    observed, without missing observations. If y contains NaNs, specifying
%    'OLS' as the method results in an error. Yule-Walker estimates are
%    available whether or not there are missing observations in y.
%
% Example:
%
%   % Create an AR(2) process from a sequence of 1000 Gaussian deviates,
%   % and assess whether the PACF is zero for lags > 2:
%
%   x = randn(1000,1);             % 1000 Gaussian deviates ~ N(0,1)
%   y = filter(1,[1 -0.6 0.08],x); % Create an AR(2) process
%   parcorr(y,'NumAR',2)           % Inspect the PACF with 95% confidence
%
% References:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
%   [2] Hamilton, J.D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also AUTOCORR, CROSSCORR, FILTER.

% Copyright 2018 The MathWorks, Inc.   

% Preprocess varargin for target axes:

try
    
    [ax,args] = internal.econ.axesparser(varargin{:});
    
catch ME
    
    throw(ME)
    
end

% This function produces a single plot:

if ~isempty(ax) && ~isscalar(ax)
    
    error(message('econ:internal:econ:axesparser:InvalidParent'));
    
end

% Parse inputs and set defaults:

if (numel(args) > 1) && isnumeric(args{2}) % Deprecated positional inputs

    % Positional input syntax:
    % [...] = parcorr(y,numLags,numAR,numSTD)

    y = args{1};
    args = args(2:end);

    try
        
        parser = checkPositionalInputs(y,args{:});
     
    catch ME
        
        throwAsCaller(ME)
     
    end
   
    numLags = parser.Results.NumLags;
    method  = parser.Results.Method;
    
else % Name-value pair inputs
    
    parser = inputParser;
    parser.addRequired ('y'      ,        @(x) validateattributes(x, {'double'}, {'vector'}, '', 'input series'));
    parser.addParameter('NumLags',    20, @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>'  0}, '', 'number of PACF lags'));
    parser.addParameter('NumAR'  ,     0, @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>=' 0}, '', 'number of AR lags'));
    parser.addParameter('NumSTD' ,     2, @(x) validateattributes(x, {'double'}, {'scalar'           '>=' 0}, '', 'number of standard deviations'));
    parser.addParameter('Method' , 'OLS', @(x) validateattributes(x, {'char' 'string'}, {'scalartext'}      , '', 'method'));
   
    y = args{1};
    args = args(2:end);
    
    try
        
        parser.parse(y,args{:});
        
    catch ME
        
        throwAsCaller(ME)
        
    end
    
    % Set lags:
 
    numLags = parser.Results.NumLags;

    if any(strcmpi('NumLags',parser.UsingDefaults)) % Default lags
       
        numLags = min(numLags,sum(~isnan(y))-1);
      
    else % User-specified lags
       
        if numLags > (sum(~isnan(y))-1)
          
            error(message('econ:parcorr:LagsTooLarge'))
        
        end
      
    end
    
    % Set method:

    if any(strcmpi('Method', parser.UsingDefaults)) % Default method
       
        if sum(~isnan(y)) < length(y)
          
            method = 'Yule-Walker';
         
        else
          
            method = 'OLS';
         
        end
      
    else % User-specified method
       
        method  = validatestring(parser.Results.Method, {'OLS' 'Yule-Walker'}, '', 'Method');
      
        if (sum(~isnan(y)) < length(y)) && strcmpi(method,'OLS')
          
            error(message('econ:parcorr:IncompatibleMissingDataMethod'))
           
        end
      
    end
   
end % Input parse

% Preprocess validated inputs:

y         = parser.Results.y;
rowSeries = (size(y,1) == 1);
y         = y(:);
N         = sum(~isnan(y)); % Effective sample size

numAR     = parser.Results.NumAR;

if numAR >= numLags
    
    error(message('econ:parcorr:NumARTooLarge'))
    
end

numSTD    = parser.Results.NumSTD;

% Compute PACF:

pacf = [1;zeros(numLags,1)];

if strcmpi(method,'OLS') % OLS estimates

% Create a lagged regression matrix & preallocate for the PACF:

    X = lagmatrix(y,1:numLags);

% Compute the PACF by fitting successive order AR models by OLS, retaining
% the last coefficient of each regression:

    for L = 1:numLags
        
        [Q,R] = qr([ones((length(y)-L),1)  X(L+1:end,1:L)],0);
        b = R\(Q'*y(L+1:end));
        pacf(L+1) = b(end);
        
   end

else % Yule-Walker estimates

    ACF = autocorr(y,numLags);

    for L = 1:numLags
        
       AR        = toeplitz(ACF(1:L))\ACF(2:(L+1));
       pacf(L+1) = AR(end);
       
    end

end

% Compute approximate confidence bounds using the approach in [1],
% equations 3.2.36 and 6.2.3, pp. 68 and 188, respectively:

bounds = [numSTD;-numSTD]./sqrt(N);
lags = (0:numLags)';

% Perform nargout-dependent operations:

nargoutchk(0,4)

if (nargout == 0) || (nargout == 4) % Create plot
    
    % Plot to gca if no parent axes is specified:

    if isempty(ax)
    
        ax = gca;
    
    end

    % Store NextPlot flag (and restore on cleanup):

    next = get(ax,'NextPlot');
    cleanupObj = onCleanup(@()set(ax,'NextPlot',next));

    % Plot the sample PACF:

    hPlot = stem(ax,lags,pacf,'filled','k-o','MarkerSize',4,'Tag','PACF');

    % Plot confidence bounds under the hypothesis that y is an AR(numAR)
    % process. The following approximation gives an indication of whether
    % the PACF is effectively zero beyond lag numAR. Confidence bounds
    % appear over the PACF only for lags greater than numAR, where the null
    % hypothesis is assumed to hold.
    
    set(ax,'NextPlot','add')
    hBounds = plot(ax,[numAR+0.5 numAR+0.5;numLags numLags],...
                      [bounds([1 1]) bounds([2 2])],'--k',...
                      'Tag','Confidence Bound');
    hXAxis = plot(ax,[0 numLags],[0 0],'-k','Tag','X Axis');

    % Return "plot object":

    h = [hPlot;hBounds;hXAxis];
    
    % Modify axes properties conditional on NextPlot flag:
    
    ax.Tag = 'PACFPlot';

    switch next
    
        case {'replace','replaceall'}
    
            %grid(ax,'on')
            %xlabel(ax,'Lag')
            ylabel(ax,'Sample Partial Autocorrelation')
            title(ax,'Sample Partial Autocorrelation Function')
            if max(pacf) <= 1
                ax.YLim(2) = 1;
            end
            
        case {'replacechildren','add'}
        
            % Do not modify axes properties
    
    end

end

if nargout > 0

    % Re-format outputs to conform to row/column orientation of y:

    if rowSeries
        
        pacf = pacf';
        lags = lags';
        bounds = bounds';
        
    end
   
end

% Suppress assignment to ans:

if (nargout > 0) && (nargout < 4)

    varargout = {pacf,lags,bounds};
    
elseif nargout == 4
    
    varargout = {pacf,lags,bounds,h};
    
end

% -------------------------------------------------------------------------
function parser = checkPositionalInputs(y,numLags,numAR,numSTD)

% Ensure the sample data is a vector:

[rows,columns] = size(y);

if (rows ~= 1) && (columns ~= 1)
    
    error(message('econ:parcorr:NonVectorInput'))
    
end

N = sum(~isnan(y)); % Effective sample size

% Ensure numLags is a positive integer or set default:

if (nargin >= 2) && ~isempty(numLags)
    
	if numel(numLags) > 1
       
        error(message('econ:parcorr:NonScalarLags'))
      
	end
   
    if (round(numLags) ~= numLags) || (numLags <= 0)
       
        error(message('econ:parcorr:NonPositiveIntegerLags'))
      
    end
   
    if numLags > (N-1)
       
        error(message('econ:parcorr:LagsTooLarge'))
      
    end
   
else
    
    numLags = min(20,N-1); % Default
   
end

% Ensure numAR is a nonnegative integer or set default:

if (nargin >= 3) && ~isempty(numAR)
    
    if numel(numAR) > 1
       
        error(message('econ:parcorr:NonScalarNumAR'))
      
    end
   
    if (round(numAR) ~= numAR) || (numAR < 0)
       
        error(message('econ:parcorr:NegativeIntegerNumAR'))
      
    end
   
    if numAR >= numLags
       
        error(message('econ:parcorr:NumARTooLarge'))
      
    end
   
else
    
    numAR = 0; % Default
   
end

% Ensure numSTD is a positive scalar or set default:

if (nargin >= 4) && ~isempty(numSTD)
    
    if numel(numSTD) > 1
       
        error(message('econ:parcorr:NonScalarSTDs'))
      
    end
   
    if numSTD < 0
       
        error(message('econ:parcorr:NegativeSTDs'))
      
    end
   
else
    
    numSTD = 2; % Default
   
end

% Set the calculation method based on the presence of missing data (NaNs).
% 'Method' is new to the name-value pair syntax; it was not available in
% the original positional syntax.

if N < length(y)
    
    parser.Results.Method  = 'Yule-Walker';
   
else
    
    parser.Results.Method  = 'OLS';
   
end

parser.Results.y       = y;
parser.Results.NumLags = numLags;
parser.Results.NumAR   = numAR;
parser.Results.NumSTD  = numSTD;