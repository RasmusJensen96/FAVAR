% Author: Fernando Ferreira
% email: fferreira@lps.ufrj.br
% Oct 2011
clear, clc

% Load data
%serie = load('./flour-price.dat')';
%serie_trval = serie(:, 1:length(serie)-30);
%serie_teste= serie(:, length(serie)-30+1:end);ee
addpath('data/sinteticos');
[serie_trval, serie_test] = sinteticas(0);

% Apply ICA (FastICA - SOBIRO)

% Load Object
TS = TimeSeries(serie_trval);

% Preprocessing
TS.preprocess();

% Create Estimator Engine
TS.assembleData();
nnobj = swirlTraining();
TS.createEstimator(nnobj);

% Apply model in test set
%load 'net_trained_TS.mat';
TS.applyModel(serie_test);


