function [Table] = GenerateSummaryStats(Data, vnames)
% Generate a table with summary statistics of data
% Rasmus M. Jensen
N = size(Data, 1);
Mean = mean(Data); Std  = std(Data); Kurt = kurtosis(Data); med = median(Data);
Min  = min(Data); Max = max(Data); Skew = skewness(Data);

Table = table(repmat(N,size(Mean, 2), 1),Mean', med', Std', Min', Max', Skew', Kurt','RowNames',vnames,'VariableNames',{'Obs','Mean','Median','Std','Minimum','Maximum','Skewness','Kurtosis'});
end
