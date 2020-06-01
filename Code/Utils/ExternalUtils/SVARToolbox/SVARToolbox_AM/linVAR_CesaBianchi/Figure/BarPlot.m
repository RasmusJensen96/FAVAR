function H = BarPlot(X, DATES)
% =======================================================================
% Creates a bar graph with positive data values stacked on the positive
% quadrant and negative data values stacked on the negative quadrant
% =======================================================================
% H = BarPlot(X)
% -----------------------------------------------------------------------
% INPUT
%   - X: data to plot [nobs x nvars]
% -----------------------------------------------------------------------
% OUTPUT
%   - H: handle to graph
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

% Modified 2020 Rasmus M. Jensen

H(1,:) = bar(DATES,(X).*(X>0),'stacked'); 
hold on;
H(2,:) = bar(DATES, (X).*(X<0),'stacked');
