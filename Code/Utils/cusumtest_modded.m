function varargout = cusumtest(varargin)
%CUSUMTEST Cusum tests for structural change
%
% Syntax:
%
%   [h,H,Stat,W,B] = cusumtest(X,y)
%   [h,H,Stat,W,B] = cusumtest(Tbl)
%   [h,H,Stat,W,B] = cusumtest(...,param,val,...)
%   [h,H,Stat,W,B,sumPlots] = cusumtest(...)
%   [h,H,Stat,W,B,sumPlots] = cusumtest(ax,...)
%   cusumtest(...)
%
% Description:
%
%   Cusum tests assess the stability of coefficients b in a multiple linear
%   regression model of the form y = X*b + e. Inference is based on a
%   sequence of sums, or sums of squares, of recursive residuals
%   (standardized one-step-ahead forecast errors) computed iteratively from
%   nested subsamples of the data. Under the null hypothesis of coefficient
%   constancy, values of the sequence outside an expected range suggest
%   structural change in the model over time.
%
% Input Arguments:
%
%   X - numObs-by-numPreds matrix of predictor data for a multiple linear                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
%       regression model.
%
%   y - numObs-by-1 vector of response data for a multiple linear
%       regression model.
%
%   Tbl - numObs-by-numPreds+1 tabular array of data for a multiple linear
%       regression model, with predictor data X in the first numPreds
%       columns and response data y in the last column.
%
%   ax - Vector of axes objects in which to plot, of length numTests. If
%       unspecified, CUSUMTEST plots to separate figures for each test.
%
%   CUSUMTEST removes observations with missing (NaN) values in the
%   predictors or the response.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'Intercept' Logical scalar or vector of length numTests indicating
%               whether or not to add an intercept when fitting the model.
%               The default is true. For any test, the number of model
%               coefficients, numCoeffs, is numPreds plus 'Intercept'.
%
%   'Test'      String vector or cell vector of character vectors, of
%               length numTests, indicating the type of test. Values are
%               'cusum' or 'cusumsq'. A 'cusum' test uses the cusum
%               statistic in [1]. A 'cusumsq' test uses the cusum of
%               squares statistic in [1]. The default is 'cusum'.
%
%   'Direction' String vector or cell vector of character vectors, of
%               length numTests, indicating the direction of iteration.
%               Values are 'forward' or 'backward'. Forward tests compute
%               recursive residuals beginning with the first numCoeffs+1
%               observations, then add observations one at a time until
%               numObs is reached. Backward tests do the same, but first
%               reverse the order of observations. The default is
%               'forward'.
%
%   'Alpha'     Scalar or vector of length numTests of nominal significance
%               levels for the tests. For cusum tests, values must be
%               between 0 and 1. For cusum of squares tests, values must
%               be between 0.01 and 0.20. The default is 0.05.
%
%   'Display'   String or character vector indicating whether or not to
%               display the results of the tests in the command window.
%               Values are 'summary' and 'off'. The default is 'off' for a
%               single test and 'summary' for multiple tests.
%
%   'Plot'      String or character vector indicating whether or not to
%               plot the results of the tests. Plots show the sequence of
%               cusums or cusums of squares, depending on the value of
%               'Test', together with critical lines determined by the
%               value of 'Alpha'. Individual plots are produced for each
%               test. Values are 'on' and 'off'. The default is 'off' when
%               CUSUMTEST is called with output arguments, 'on' otherwise.
%
%   The number of tests, numTests, is determined by the length of any
%   vector parameter value. Scalar or single string parameter values are
%   expanded to vectors of length numTests. Vector values must all have
%   length numTests. If any parameter value is a row vector, so is output
%   h (array outputs retain their specified dimensions).
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to
%       numTests. Hypotheses are independent of the value of 'test':
%
%       H0: Coefficients in b are equal in all sequential subsamples.
%
%       H1: Coefficients in b have changed during the period of the sample.
%
%       Values of h equal to true indicate rejection H0 in favor of H1.
%       Values of h equal to false indicate a failure to reject H0.
%
%   H - Array of Boolean decisions for the sequence of tests, of size equal
%       to numTests-by-(numObs-numPreds). For forward tests, column indices
%       correspond to times numPreds+1, ..., numObs. For backward tests,
%       column indices correspond to times numObs-(numPreds+1), ..., 1.
%       Rows corresponding to tests where 'Intercept' is true contain one
%       less iteration, and the value in the first column defaults to
%       false. The output h is any(H,2).
%
%   Stat - Array of statistics, the same size as H. Values in any row
%       depend on the value of 'Test' (cusum or cusum of squares).
%       Indexing corresponds to the indexing in H. Rows corresponding to
%       tests where 'Intercept' is true contain one less iteration, and the
%       value in the first column defaults to NaN.
%
%   W - Array of standardized recursive residuals, as in [1], the same size
%       as H. Indexing corresponds to the indexing in H. Rows corresponding
%       to tests where 'Intercept' is true contain one less iteration, and
%       the value in the first column defaults to NaN.
%
%   B - Array of recursive coefficient estimates, of size (numPreds+1)-by-
%       (numObs-numPreds)-by-numTests. Index (i,j,k) corresponds to
%       coefficient i at iteration j for test k. At iteration j, the test
%       computes the coefficients 
%
%       b = X(1:numPreds+j,inRegression)\y(1:numPreds+j),
%
%       where inRegression is a logical vector indicating the predictors in
%       the regression at iteration j.
%
%       Initially constant predictors are held out of forward regressions
%       until their data changes, if they introduce multicollinearity (see
%       Notes). Coefficient estimates for hold-out iterations default to
%       NaN. Similarly, terminally constant predictors are held out of
%       backward regressions.
%
%       Tests where 'Intercept' is true contain one less iteration, and the
%       value in the first column defaults to NaN. Tests where 'Intercept'
%       is false contain one less coefficient, and the value in the first
%       row, corresponding to the intercept, defaults to NaN.
%
%   sumPlots - Array of handles to plotted graphics objects, of size
%       3-by-numTests. This output is available only if the 'Plot'
%       parameter is 'on'.
%
% Notes:
%
%   o Cusum tests provide useful diagnostics for a variety of model
%     misspecifications, including gradual structural change, multiple
%     structural changes, missing predictors, and neglected nonlinearities.
%     The tests, however, have little power to detect structural changes
%     late in the sample period, or multiple changes that produce
%     cancellation in the cusums.
%
%   o The cusum of squares test is a "useful complement to the cusum
%     test, particularly when the departure from constancy of the
%     [recursive coefficients] is haphazard rather than systematic" [1]. It
%     has greater power in cases where multiple shifts may cancel out. It
%     is often suggested for detecting structural breaks in volatility.
%
%   o CUSUMTEST handles initially constant predictor data using the method
%     suggested in [1]. If a predictor's data is constant for the first
%     numCoeffs observations and this results in multicollinearity with an
%     intercept or another predictor, then CUSUMTEST drops the predictor
%     from regressions and the computation of recursive residuals until its
%     data changes. Similarly, CUSUMTEST temporarily holds out terminally
%     constant predictors from backward regressions. Initially constant
%     predictors in backward regressions, or terminally constant predictors
%     in forward regressions, are not held out by CUSUMTEST, and can lead
%     to rank deficiency in terminal iterations.
%
%   o Critical lines used for inference are computed in essentially
%     different ways for the two test statistics. For cusums, the normal
%     cdf equation in [1] is solved dynamically for each value of 'Alpha'.
%     For cusums of squares, parameter values are interpolated from the
%     table in [2], using the method suggested in [1]. Sample sizes with
%     degrees of freedom less than 4 are below tabulated values, and
%     critical lines cannot be computed. Sample sizes with degrees of
%     freedom greater than 202 are above tabulated values, and the critical
%     value associated with the largest tabulated sample size is used.
%
%   o Nominal significance levels are specified by 'Alpha'. The actual size
%     of a test depends on a number of assumptions and approximations used
%     to compute the critical lines. Plots of the recursive residuals are
%     the best indicator of structural change, and [1] suggests that the
%     tests "should be regarded as yardsticks for the interpretation of
%     data rather than leading to hard and fast decisions."
%
%   o Basic diagnostic plots of the recursive coefficient estimates are
%     produced by plot(B(:,:,n)') for test n. Similar plots, with robust
%     standard error bands, are produced by RECREG.
%
% Example:
%
%   % Test the stability of an explanatory model of real GNP for a period
%   % spanning World War II:
%
%   load Data_NelsonPlosser
%   span = (1915 <= dates) & (dates <= 1970);
%   Mdl = DataTable(span,[4,5,10,1]);
%   h = cusumtest(Mdl,'Plot','on')
%
% References:
%
%   [1] Brown, R. L., J. Durbin, and J. M. Evans. "Techniques for Testing
%       the Constancy of Regression Relationships Over Time." Journal of
%       the Royal Statistical Society. Series B, Vol. 37, 1975, pp.
%       149-192.
%
%   [2] Durbin, J. "Tests for Serial Correlation in Regression Analysis
%       Based on the Periodogram of Least Squares Residuals." Biometrika.
%       Vol. 56, 1969, pp. 1-15.
%
% See also CHOWTEST, RECREG, FITLM

% Copyright 2015 The MathWorks, Inc.

% Preprocess varargin for target axes.

try
    
    [targetAxes,args] = internal.econ.axesparser(varargin{:});
    
catch ME
    
    throw(ME)
    
end

areTargets = ~isempty(targetAxes);

% Parse inputs and set defaults.

Data = args{1};

% Handle dataset array inputs:

if isa(Data,'dataset')  
    try
        Data = dataset2table(Data);
    catch 
        error(message('econ:cusumtest:DataNotConvertible'))
    end
end

parseObj = inputParser;

addRequired(parseObj,'Data',...
            @(x)validateattributes(x,{'numeric','table','dataset'},...
                                     {'2d','nonempty'}))

if isnumeric(Data)
    
    addRequired(parseObj,'y',...
        @(x)validateattributes(x,{'numeric'},{'vector','nonempty'}))
    
elseif size(Data,2) < 2 % Table without response
    
    error(message('econ:cusumtest:DataWithoutResponse'))
    
end   

addParameter(parseObj,'Intercept',true,...
        @(x)validateattributes(x,{'numeric','logical'},{'vector','binary'}))
    
addParameter(parseObj,'Test','cusum',...
        @(x)validateattributes(x,{'char','cell','string'},{'vector'}))
    
addParameter(parseObj,'Direction','forward',...
        @(x)validateattributes(x,{'char','cell','string'},{'vector'}))

addParameter(parseObj,'Alpha',0.05,...
        @(x)validateattributes(x,{'numeric'},{'vector','>',0,'<',1}))
    
addParameter(parseObj,'Display','off',...
        @(x)validateattributes(x,{'char','string'},{'vector'}))

addParameter(parseObj,'Plot','off',...
        @(x)validateattributes(x,{'char','string'},{'vector'}))

parse(parseObj,args{:});

if isnumeric(Data)
    
    y = parseObj.Results.y;
    
end

iFlag = logical(parseObj.Results.Intercept);

whichTest = parseObj.Results.Test;
if iscell(whichTest) || (isstring(whichTest) && ~isscalar(whichTest))
    for testIdx = 1:length(whichTest)
        whichTest{testIdx} = ...
            validatestring(whichTest{testIdx},{'cusum','cusumsq'});
    end
    whichTest = cellstr(whichTest);
else
    whichTest = validatestring(whichTest,{'cusum','cusumsq'});
end

whichDirection = parseObj.Results.Direction;
if iscell(whichDirection) || (isstring(whichDirection) && ~isscalar(whichDirection))
    for dirIdx = 1:length(whichDirection)
        whichDirection{dirIdx} = ...
            validatestring(whichDirection{dirIdx},{'forward','backward'});
    end
    whichDirection = cellstr(whichDirection);
else
    whichDirection = validatestring(whichDirection,{'forward','backward'});
end

alpha = parseObj.Results.Alpha;

dispFlag = validatestring(parseObj.Results.Display,{'off','summary'});

plotFlag = validatestring(parseObj.Results.Plot,{'off','on'});

% Check parameter dimensions, expand scalars and single strings:

[numTests,rowOutput,iFlag,whichTest,whichDirection,alpha] = ...
    sizeCheck(iFlag,whichTest,whichDirection,alpha);

% Load critical values for cusumsq test, as necessary:

needDurbinData = any(strcmp(whichTest,'cusumsq'));
if needDurbinData

    CV = load('Data_Durbin');
    C = CV.C; % Table 1 from [2].
    
end

% Set 'Display' parameter:

if any(strcmp(parseObj.UsingDefaults,'Display')) 
    if numTests == 1
        dispFlag = 'off';
    else
        dispFlag = 'summary';
    end
end

needDisplay = strcmp(dispFlag,'summary');
if needDisplay
    
    fprintf('\nRESULTS SUMMARY\n\n')
    
end

% Set 'Plot' parameter:

if any(strcmp(parseObj.UsingDefaults,'Plot'))
    if nargout > 0
        plotFlag = 'off';
    else
        plotFlag = 'on';
    end
end

needPlot = strcmp(plotFlag,'on');

% Check target axes for appropriate size:

if needPlot && areTargets && ...
   (~isvector(targetAxes) || (length(targetAxes) ~= numTests))
    
        error(message('econ:cusumtest:TargetsWrongSize'))
    
end

% Get X, y:

if isnumeric(Data)
    
    X = Data;
    
    if length(y) ~= size(X,1)
             
    	error(message('econ:cusumtest:ResponseVectorWrongSize'))
          
    else 
    
        y = y(:);
        
    end
    
else
    
    X = table2array(Data(:,1:(end-1)));
    y = table2array(Data(:,end));
    
end

% Remove observations with missing values:

D = [X,y];
D(any(isnan(D),2),:) = [];
X = D(:,1:end-1);
y = D(:,end);
[numObs,numPreds] = size(X);

% Check for sufficient observations:

dfe0 = numObs-numPreds;
isIntercept = any(iFlag);

if dfe0 < isIntercept
        
	error(message('econ:cusumtest:TooFewObservations'))
         
end

% Preallocate outputs:

h = false(numTests,1);
H = false(numTests,dfe0);
Stat = NaN(numTests,dfe0);
W = NaN(numTests,dfe0);
B = NaN(numPreds+1,dfe0,numTests);

if needPlot
    
    sumPlots = gobjects(3,numTests);
    
end

% Perform the tests:

for t = 1:numTests
    
    testIntercept = iFlag(t);
    testType = whichTest{t};
    testDirection = whichDirection{t};
    testAlpha = alpha(t);
    
    % Prepare data:
    
    if strcmp(testDirection,'forward')
        
        XTest = X;
        yTest = y;
        
    else
        
        XTest = flipud(X);
        yTest = flipud(y);
        
    end
    
    inRegression = true(size(X));
    numCoeffs = numPreds;
    dfe = dfe0;
 
    for j = 1:numPreds

        firstChange = find(XTest(:,j) ~= XTest(1,j),1);
        
        if isempty(firstChange)
            
            inRegression(:,j) = false;
            
        else
        
            inRegression(1:firstChange-1,j) = false;
        
        end

    end

    if testIntercept

        XTest = [ones(numObs,1),XTest]; %#ok
        
        inRegression = [true(numObs,1),inRegression]; %#okb
        numCoeffs = numCoeffs+1;
        dfe = dfe0-1;
        
    else
        
        % Allow constant predictors if no multicollinearity:
        
        for i = 1:numObs
            
            if sum(inRegression(i,:) == false) < 2
                
                inRegression(i,:) = true;
                
            end
            
        end

    end

    % Run the regressions:
        
    X0 = XTest(1:numCoeffs,inRegression(numCoeffs,:));        
    y0 = yTest(1:numCoeffs);
    b0 = X0\y0;
    res0 = y0-X0*b0;
    SSE0 = res0'*res0;
    
	for r = numCoeffs+1:numObs
        
        jIter = r-numPreds; % Write to this column
        
        % Compute recursive residuals:
        
        x1 = XTest(r,inRegression(r-1,:)); 
        y1 = yTest(r);        
        w1 = (y1-x1*b0)./sqrt(1+(x1/(X0'*X0))*x1');
        W(t,jIter) = w1;
        
        % Compute coefficient estimates:
        
        X0 = XTest(1:r,inRegression(r,:));        
        y0 = [y0;y1]; %#ok
        b0 = X0\y0;
        
        if testIntercept
            
            bRows = inRegression(r,:); 
            
        else
            
            bRows = [false,inRegression(r,:)];
            
        end
        
        B(bRows,jIter,t) = b0;

	end % End regression loop
    
    % Compute statistics:
    
    j1 = testIntercept+1; % First column with stats
    
    SSE = cumsum([SSE0,W(t,j1:end).^2]);

    if strcmp(testType,'cusum')
        
        RMSE = sqrt(SSE(end)/dfe);
        Stat(t,j1:end) = cumsum(W(t,j1:end))/RMSE;

    else % Cusum of squares

        Stat(t,j1:end) = SSE(2:end)/SSE(end);

    end

    % Compute inferences:
    
    s = Stat(t,j1:end);
    iter = 1:dfe;
    
	if strcmp(testType,'cusum')
        
        % Compute parameter a:
        
        Q = @(a)(normcdf(-a));
        paramEqn = @(a)(Q(3*a) + exp(-4*a^2)*(1-Q(a)) - testAlpha/2);
        options = optimoptions('fsolve','Display','off');
        a = fsolve(paramEqn,0,options);
        
        % Compute critical lines:
        
        rdfe = sqrt(dfe);
        c1 = (2*a/rdfe)*iter + a*rdfe;
        c2 = -(2*a/rdfe)*iter - a*rdfe;
        H(t,j1:end) = (s > c1) | (s < c2);
       
	else % Cusum of squares
        
        % Compute parameter c0:
        
        c0 = DurbinTableLookup(testAlpha,dfe,C);
        
        % Compute critical lines:
        
        c1 = c0 + iter/dfe;
        c2 = -c0 + iter/dfe;
        
	end
    
    % Perform the test:
    
    H(t,j1:end) = (s > c1) | (s < c2);    
    testDecision = any(H(t,:));
    h(t) = testDecision;
    
    % Update display:

    if needDisplay
        
        if testIntercept

            testInterceptString = 'yes';

        else

            testInterceptString = 'no';

        end

        if testDecision

            testDecisionString = 'Reject coefficient stability';

        else

            testDecisionString = 'Fail to reject coefficient stability';

        end

        fprintf('***************')
        fprintf('\nTest %d\n',t)

        fprintf('\nTest type: %s',testType)
        fprintf('\nTest direction: %s',testDirection)
        fprintf('\nIntercept: %s',testInterceptString)
        fprintf('\nNumber of iterations: %d\n',dfe)

        fprintf('\nDecision: %s',testDecisionString)
        fprintf('\nSignificance level: %.4f\n\n',testAlpha)

    end

    % Create plots:

    if needPlot
       
        if dfe <= 1
            
            warning(message('econ:cusumtest:NoPlot'))
            
        else
            
            % Set axes:
            
            if areTargets
                
                ax = targetAxes(t);
                
            else
                
                figure
                ax = gca;
                
            end
            
            % Store NextPlot flag (and restore on cleanup):

            next = get(ax,'NextPlot');
            cleanupObj = onCleanup(@()set(ax,'NextPlot',next));
            
            % Create plot:
            
            hPlot = plot(ax,s,'b-','LineWidth',2,'Tag','CusumPlot');

            set(ax,'NextPlot','add');
            hUpper = plot(ax,c1,'r--','Tag','UpperCriticalLine');
            hLower = plot(ax,c2,'r--','Tag','LowerCriticalLine');
            
            sumPlots(:,t) = [hPlot;hUpper;hLower];
            
            % Modify other axes properties conditional on NextPlot flag:

            ax.Tag = ['CusumPlot ',num2str(t)];

            switch next

                case {'replace','replaceall'}

                    xlabel(ax,'Iteration')
                    ylabel(ax,testType)
                    title(ax,['{\bf Test ',num2str(t),'}'])
                    legend(ax,[hPlot,hUpper],{'Statistic',...
                           'Critical Lines'},'Location','best');
                    axis(ax,'tight')
                    grid(ax,'on')

                case {'replacechildren','add'}

                    % Do not modify axes properties

            end
            
        end
        
    end
    
end % End testing loop

% Display vector row output:

if rowOutput
    
    h = h';
    
end

% Suppress assignment to ans:

nargoutchk(0,8)

if (nargout > 0) && (nargout < 8)
    
    varargout = {h,H,Stat,W,B,c1,c2};

elseif (nargout == 8) && needPlot
    
    varargout = {h,H,Stat,W,B,c1,c2,sumPlots};
    
elseif (nargout == 8) && ~needPlot
    
    error(message('econ:cusumtest:NoPlotHandles'))
    
end

%-------------------------------------------------------------------------
% Check parameter dimensions, expand scalars and single strings
function [numTests,rowOutput,varargout] = sizeCheck(varargin)

numCheck = length(varargin);
varargout = cell(1,numCheck);

checkLengths = zeros(numCheck,1);

% Find parameter lengths:
for i = 1:numCheck
    
    ivar = varargin{i};
    
    if isnumeric(ivar) || islogical(ivar) || iscell(ivar)
        
        checkLengths(i) = length(ivar);
        
    else
        
        checkLengths(i) = 1; % Single string
        varargin{i} = varargin(i); % Convert to cell
        
    end
    
end

% Check nonscalar parameter lengths for commensurability, set numTests:

nonScalarIdx = find(checkLengths > 1);

if ~isempty(nonScalarIdx)
    
    ns1 = nonScalarIdx(1);
    
    if any(checkLengths(nonScalarIdx) ~= checkLengths(ns1))
        
        error(message('econ:cusumtest:ParameterSizeMismatch'))

    end
    
    numTests = checkLengths(ns1);
    
else
    
    numTests = 1;

end

% Set row output:

[~,n] = cellfun(@size,varargin);
rowOutput = any(n > 1);

% Assign output:

for i = 1:numCheck

    if checkLengths(i) == 1

        varargout{i} = repmat(varargin{i},numTests,1);

    else

        varargout{i} = varargin{i};

    end

end

%-------------------------------------------------------------------------
% Get critical value from Durbin's table
function c0 = DurbinTableLookup(testAlpha,dfe,C)

levels = [0.10 0.05 0.025 0.01 0.005];
degrees = [(1:60) (62:2:100)];

a = testAlpha/2;
n = (dfe/2)-1;

if (a > levels(1)) || (a < levels(end))

    error(message('econ:cusumtest:AlphaOutOfRange'))

elseif n < degrees(1)

    error(message('econ:cusumtest:SmallDfe'))

elseif n > degrees(end)

    degrees(end) = n; % Treat sample size as asymptotic

end

c0 = interp2(levels,degrees,C,a,n,'linear');