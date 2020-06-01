function [Ytph_Filtered, FE, ytph] = RandomwalkForecast(Y, t, h)
% Function returning forecast error and filtered value of a random walk with
% no drift.
% Rasmus M. Jensen 2020

Ytph_Filtered = Y(t-1);
ytph = Y(t-1+h);
FE = ytph-Ytph_Filtered;
end

