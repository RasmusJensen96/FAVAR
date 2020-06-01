% FAVAR identification, assesment and forecasting 
% Rasmus M. Jensen 2020
clear all;
tic;
addpath(genpath('Data'));
addpath(genpath('Utils'));
addpath(genpath('SVARToolbox')); clear all; clc;
set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
signRest         = 0;  % Aux switch indicating whether or not sign-restrictions should be applied
FFR1d            = 0;  % Switch to determine whether FFR should be assumed I(1) or levels I(0)
Pre2009          = 0;  % Aux switch, 0 = 1959:2019, 1 = 1959:2008
type             = 2;  % Switch 1 = PCA extraction of factors, 2 = Kalman Smoothing of factors. NB. Kalman is computationally very heavy.
nsteps           = 24; % Number of steps IRF
ConfLevel        = 68; % Confidence level percentage for bootstrap
Information_Crit = 0; % Switch 1: Lag length and numfacs by IC
var_numbers = [1, 19, 24, 66, 73, 80, 128, 129, 130 131];
Order = [1 2 3];
%% Load 1983:2019 Data
%-----------------------------Load Data------------------------------------
% PrepData is a seperate script importing and preparing all dataseries,
% including removing outliers, standardizing and running the EM-algorithm
% it outputs all data neccesary for estimating the FAVAR, as well as
% an auxiliary vector of sample dates.
PrepData
%%
plot(ydata(:,3));
%%
t = size(ydata, 1);
[~,~,r,~] = regress(ydata(:,3),[ones(t,1)'; 1:t]')
%%
plot(r);
hold on;
plot(ydata(:,3))