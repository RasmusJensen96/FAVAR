% Script collecting data and calculating the real (deflated) FFR
url = 'https://fred.stlouisfed.org/';
c = fred(url);
series = 'PCEPILFE';
d = fetch(c,series);
y = d.Data(:,2);
date = datetime(d.Data(:,1),'ConvertFrom','datenum');
for i = 13:size(y, 1)
    y_c(i,:) = (y(i)-y(i-12)); % YoY-PCE
end
y_c = y_c(13:end,:);
date = date(13:end);
series = 'FEDFUNDS';
d = fetch(c,series);
ffr = d.Data(1:end-1,:); 
ffr = ffr(55+12:end,:);
rffr = ffr(:,2) - y_c;
plot(date,rffr); hold on; plot(date, ffr(:,2));
yline(0,':')
%% Calculating the real federal funds rate:
clearvars -except rffr date
save("rffr")
