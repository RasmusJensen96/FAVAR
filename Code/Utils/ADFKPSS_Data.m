function tests = ADFKPSS_Data(Y, namesY, lags)
for i=1:size(Y, 2)
    [~,pValue,stat,cValue,reg] = kpsstest(Y(:,i),'lags',1:12);
    eval(join(['tests.kpss.', char(namesY{i}),' = [stat;cValue;pValue; reg.DWStat]']));
    [~,pValue,stat,cValue,reg] = adftest(Y(:,i),'lags',1:12);
    eval(join(['tests.adf.', char(namesY{i}),' = [stat;cValue;pValue; reg.DWStat]']));
end
end