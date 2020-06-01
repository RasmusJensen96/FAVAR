function [SeasAdj] = Deseasonalize12StableFilter(vY)
% Seasonally adjusts an input series using a 12 month stable seasonal filter
% Input = vY: Vector of unadjusted time-series
% Rasmus M. Jensen 2020
T = length(vY);
%%
sW13 = [1/24;repmat(1/12,11,1);1/24];
yS = conv(IP_NA,sW13,'same'); % Note what we require the signal proc. TB
yS(1:6) = yS(7); yS(T-5:T) = yS(T-6); % To not abandon any observatitons we let the MA trend oof the extremum observations equal the first filtered value
xt = IP_NA-yS; %% The stable seasonal component
%% 
s = 12;
sidx = cell(s,1);
for i = 1:s
 sidx{i,1} = i:s:T;
end
%% 
sst = cellfun(@(x) mean(xt(x)),sidx);
% Put smoothed values back into a vector of length N
nc = floor(T/s); % no. complete years
rm = mod(T,s); % no. extra months
sst = [repmat(sst,nc,1);sst(1:rm)];
% Center the seasonal estimate (additive)
sBar = mean(sst); % for centering
sst = sst-sBar;
%% 
SeasAdj = IP_NA - sst; 
end