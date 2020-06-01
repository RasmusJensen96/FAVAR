function [Ytph_Filtered FE] = MeanModelFC(Y, t,WinSize, h);

yt = Y(t-WinSize-1:t-1);
ytph = Y(t+h-1);

Ytph_Filtered = mean(yt);

FE = ytph - Ytph_Filtered;