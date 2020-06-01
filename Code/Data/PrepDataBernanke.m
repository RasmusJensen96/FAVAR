% Script loading and preparing all data:
% Full
%%
% Factor dataset:
xdata = xlsread("FactorData.xlsx");
% load data on inflation, unemployment and interest rate 
ydata = readmatrix("VAR_Data.xlsx"); % Y_data transformed in sep. file:
% load transformation codes (see file transx.m)
tcode = readmatrix("Transform_codes.xlsx");
tcode = [tcode,5,5,2];
slowcode = readmatrix("Slow_Code.xlsx");
% load the file with the dates of the data (Monthly)
yyearlab = readmatrix("DatesY.xlsx");
ydate = datetime(yyearlab,'ConvertFrom','datenum');
ydate = ydate + 365.25 * 1900;
xyearlab = readmatrix("DatesX.xlsx");
xdate = datetime(xyearlab,'ConvertFrom','datenum');
xdate = xdate + 365.25 * 1900;
samplendY    = ydate <= datetime("16-Aug-2001");
samplendX    = xdate <= datetime("16-Aug-2001");
DATESnm = yyearlab(samplendY);
DATES = ydate(samplendY);

%xindex = [find(contains(namesXY,'FEDFUNDS')); find(contains(namesXY,'INDPRO')); find(contains(namesXY,'CPIAUCSL'))];
xdata = xdata(samplendX,:);
ydata = ydata(samplendY,1:3); % Inflation, Industrial production, FFR
ydata = ydata(:,Order);
ydata_st = (ydata - mean(ydata))./std(ydata);
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
xt=prepare_missing(xdata,tcode);
% Reduce sample to usable dates due to transformations.
xt=xt(3:end,:); 

% Remove outliers using auxiliary function remove_outliers() (From: McCracken 2015)
% data = matrix of transformed series with outliers removed
% n = number of outliers removed from each series
[data,n]=remove_outliers(xt);
%%
% jj = Bai & Ng Information criterios to use 1:3
jj=2;
% kmax = maximum number of factors to include in the EM-algorithm
kmax=10;
% Demean & Standardize.
DEMEAN = 2;
[~,~,~,~,x2] = factors_em(data,kmax,jj,DEMEAN); % EM-Algorithm for NA's  (From: McCracken 2015)
x2 = x2(:,xindex(1:end-3));
% Demeaned Dataset
x2_st = (x2-mean(x2))./std(x2);
slowcode = slowcode(xindex(1:end-3));
tcode = tcode(xindex);
%%
clearvars jj kmax DEMEAN n data xt samplestartX samplestartY samplendY ...
     yearlab ydate xdate yyearlab xyearlab datefinx xindex yindex samplendX...
     xdata namesX