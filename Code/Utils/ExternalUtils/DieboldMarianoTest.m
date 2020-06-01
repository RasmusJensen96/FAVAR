function [dmteststat] = DieboldMarianoTest(FE1, FE2);
% Infer equality of two vectors of forecasts qua the Diebold Mariano 1995
% test, one-period ahead.
% Rasmus M. Jensen 2020
Loss = FE1 - FE2;
mLoss = mean(Loss);
vLoss = var(Loss);
dmteststat = mLoss * sqrt(1/size(FE1, 1) * vLoss)^(-1);
