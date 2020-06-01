% Script loading and preparing all data with a split  post-Volcker:
% Factor dataset:
xdata = xlsread("FactorData.xlsx");
% load data on inflation, unemployment and interest rate 
ydata = readmatrix("VAR_Data.xlsx"); % Y_data transformed in sep. file:
% load transformation codes (see file transx.m)
tcode = readmatrix("Transform_codes.xlsx");
slowcode = readmatrix("Slow_Code.xlsx");
% load the file with the dates of the data (Monthly)
yyearlab = readmatrix("DatesY.xlsx");
ydate = datetime(yyearlab,'ConvertFrom','datenum');
ydate = ydate + 365.25 * 1900;
xyearlab = readmatrix("DatesX.xlsx");
xdate = datetime(xyearlab,'ConvertFrom','datenum');
xdate = xdate + 365.25 * 1900;
samplendY    = ydate <= max(xdate);
DATESnm = yyearlab(samplendY);
DATES = ydate(samplendY);
ydata = ydata(samplendY,1:3); % Inflation, Industrial production, FFR
if FFR1d == 1
    ydata(:,3) = [NaN; ydata(2:end,3)-ydata(1:end-1,3)];
end
ydata = ydata(:,Order);
ydata = ydata(3:end,:);
DATES = DATES(3:end,:);

namesX = readcell("FactorNames.xlsx")';
namesY = [{'CPI'}; {'IP'} ; {'FFR'}];
namesY = namesY(Order);
namesXY = [namesX ; namesY];
yindex  = logical(sum(string(namesXY) == ["FEDFUNDS","INDPRO","CPIAUCSL"],2));
xindex  = logical(ones(size(yindex,1),1) - yindex);
xindex([58, 69]) = 0; %% only NA-s in the pre-sample
namesXY = namesXY(xindex);
xdata = xdata(:,xindex(1:end-3));
slowcode = slowcode(:,xindex(1:end-3));
%% Two subsamples
YVolcker = DATES <= datetime("31-Dec-1983");
XVolcker = xdate <= datetime("31-Dec-1983");
PreDates = DATES(YVolcker);
ydataPre = ydata(YVolcker,:);
xdataPre = xdata(XVolcker,:);
PreDatesNum = DATESnm(YVolcker);
YVolcker = DATES > datetime("31-Oct-1983");
XVolcker = xdate > datetime("31-Oct-1983");
ydataPost = ydata(YVolcker,:);
ydataPost = ydataPost(3:end,:);
xdataPost = xdata(XVolcker,:);
PostDates = DATES(YVolcker);
PostDatesNum = DATESnm(YVolcker);
ydataPre_st = (ydataPre - mean(ydataPre))./std(ydataPre);
ydataPost_st = (ydataPost - mean(ydataPost))./std(ydataPost);
namesY = [{'CPI'} {'IP'} {'FFR'}];
namesY = namesY(Order)';
%%
xtPre=prepare_missing(xdataPre,tcode);
xtPost=prepare_missing(xdataPost,tcode);
% Reduce sample to usable dates due to transformations.
xtPre=xtPre(3:end,:);
xtPost=xtPost(3:end,:); 

% Remove outliers using auxiliary function remove_outliers() (From: McCracken 2015)
% data = matrix of transformed series with outliers removed
% n = number of outliers removed from each series
[dataPre,~]=remove_outliers(xtPre);
[dataPost,~]=remove_outliers(xtPost);
%%
% jj = Bai & Ng Information criterios to use 1:3
jj=2;
% kmax = maximum number of factors to include in the EM-algorithm
kmax=10;
% DEMEAN = 0, no demean, 1, demean, 2 demean & standardize, 3 Recurs.
% Demean & Standardize.
DEMEAN = 2;
%%
[~,~,~,~,x2Pre] = factors_em(dataPre,kmax,jj,DEMEAN); % EM-Algorithm for NA's  (From: McCracken 2015)
%%
[~,~,~,~,x2Post] = factors_em(dataPost,kmax,jj,DEMEAN); % EM-Algorithm for NA's  (From: McCracken 2015)
%%
% Demeaned Dataset
x2Pre_st = (x2Pre-mean(x2Pre))./std(x2Pre);
x2Post_st = (x2Post-mean(x2Post))./std(x2Post);
%%
clearvars jj kmax DEMEAN dataPost dataPre tcode DATES yyearlab YVolcker...
    ydate xyearlab XVolcker xtPre xtPost xdate xdata ydata DATESnm sampleendY...
    xdataPost xdataPre
    