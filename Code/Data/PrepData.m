% Script loading and preparing all data:
% Sample
%%
% Factor dataset:
xdata = xlsread("FactorData.xlsx");
% load data on inflation, unemployment and interest rate 
ydata = readmatrix("VAR_Data.xlsx"); % Y_data transformed in sep. file:
%ydata(:,3) = [NaN; ydata(2:end,3)-ydata(1:end-1,3)]; % 1. diff of fed_funds
% load transformation codes (see file transx.m)
tcodes = readmatrix("Transform_codes.xlsx");
slowcode = readmatrix("Slow_Code.xlsx");
% load the file with the dates of the data (Monthly)
yyearlab = readmatrix("DatesY.xlsx");
ydate = datetime(yyearlab,'ConvertFrom','datenum');
ydate = ydate + 365.25 * 1900;
xyearlab = readmatrix("DatesX.xlsx");
xdate = datetime(xyearlab,'ConvertFrom','datenum');
xdate = xdate + 365.25 * 1900;
samplestartX  = xdate >= '01-Nov-1983';
samplestartY  = ydate >= '01-Jan-1984';
samplendY    = ydate <= max(xdate);
DATESnm = yyearlab(samplestartY & samplendY);
DATES = ydate(samplestartY & samplendY); %% Sample date-vector
%%
if FFR1d == 1;
    ydata(:,3) = [NaN; ydata(2:end,3)-ydata(1:end-1,3)];
    ydata = ydata(2:end,:); %DATES = DATES(2:end,:);DATESnm = DATESnm(2:end,:);
    samplestartX = samplestartX(2:end);
end

ydata = ydata(samplestartY & samplendY,1:3); % Inflation, Industrial production, FFR
ydata = ydata(:,Order);
% t = size(ydata, 1);
% [~,~,r,~] = regress(ydata(:,3),[ones(t,1)'; 1:t]'); Blanchard-Quah style
% detrending 
% r_std = r./std(r);
%ydata = ydata(:,Order);
ydata_st = (ydata - mean(ydata))./std(ydata);
% ydata(:,3) = r; 
% ydata_st(:,3) = r_std;
xdata = xdata(samplestartX,:);
%%
% load the file with the names of the variables. 
namesX = readcell("FactorNames.xlsx")';
namesY = [{'CPI'}; {'IP'} ; {'FFR'}];
namesY = namesY(Order);
namesXY = [namesX ; namesY];
yindex  = logical(sum(string(namesXY) == ["FEDFUNDS","INDPRO","CPIAUCSL"],2));
xindex  = logical(ones(size(yindex,1),1) - yindex);
namesXY = namesXY(xindex);
%%
xt=prepare_missing(xdata,tcodes);
% Reduce sample to usable dates due to transformations.
xt=xt(3:end,:); %% Start from january 1983

% Remove outliers using auxiliary function remove_outliers() (From: McCracken 2015)
% data = matrix of transformed series with outliers removed
% n = number of outliers removed from each series
[data,n]=remove_outliers(xt);
%%
% jj = Bai & Ng Information criteria to use 1:3
jj=2;
% kmax = maximum number of factors to include in the EM-algorithm
kmax=10;
% DEMEAN = 0, no demean, 1, demean, 2 demean & standardize, 3 Recurs.
% Demean & Standardize.
DEMEAN = 2;
[~,~,~,~,x2] = factors_em(data,kmax,jj,DEMEAN); % EM-Algorithm for NA's  (From: McCracken 2015)
x2 = x2(:,xindex(1:end-3));
% Demeaned Dataset
tcodes = tcodes(xindex(1:end-3));
slowcode = slowcode(xindex(1:end-3));
x2_st = (x2-mean(x2))./std(x2);
%%
clearvars jj kmax DEMEAN n data xt samplestartX samplestartY samplendY ...
    yearlab ydate xdate yyearlab xyearlab datefinx r r_std t yyearlab ...
    yindex xindex xdata namesX