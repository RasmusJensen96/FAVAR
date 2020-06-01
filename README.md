# Structural FAVAR
This repository hosts a framework for implementing a factor-augmented vector-autoregression (FAVAR)-model in MatLab. In addition it contains replication codes for my personal Master Thesis.

## Structure
For replicating figures and tables from the thesis, several scripts are utilized:
* 'Main_PostVolcker.m'  Constructs the structural FAVAR for the subsample 01:1984-08:2019
* 'Main_Bernanke.m'     Constructs the FAVAR for the subsample 01:1959-08:2001. 
* 'Main_Forecast.m'     Uses the reduced-form FAVAR to conduct _pseudo_-Out-of-Sample forecasts and hypothesis tests
* 'IRFplotsscript.m'    Plots the IRFs estimated from the Main-scripts. This script requires two additional scripts to be executed:
  - 'Main_Fullsample.m' Estimates IRFs from the full 01:1959-08:2019-sample
  - 'Main_Threshold.m'  Generates IRFs from the threshold FAVAR on the samples 01:1959-12:1983 and 01:1984-08:2019

Several options are available in each script.
To recreate the Bernanke et. al (2005)-paper and section 2 of the thesis:
* Following options should be set in 'Main_Bernanke.m':
  - plag = 13;      % Number of VAR-lags
  - faclag = 1      % Number of lags in the measurement equation
  - type = 1        % PCA-extraction of factors (1 = Bernanke, 2 = section 2 of the thesis)
  - K = 3           % Number of primitive factors
  - FFR1d = 0       % Consider the federal funds rate I(0) 
  - Cholesky = 1    % Recursive identification
  - Order = [1,2,3] % Infl, IP, FFR
  
To recreate the main model used in the Thesis, section 3.
* Following options should be set in 'Main_PostVolcker.m':
  - plag = 3;       % Number of VAR-lags
  - faclag = 3      % Number of lags in the measurement equation
  - type = 2        % Kalman Smooother-extraction of factors
  - K = 3           % Number of primitive factors
  - FFR1d = 1       % Consider the federal funds rate I(1) 
  - Cholesky = 0    % Blanchard Quah (1989) identification
  - Order = [2,3,1] % IP, FFR, Infl
  
  ## Tables and figures in the thesis
  
* 'StationarityScript.m' - generates the following tables/tests
  - Johansen-test table 
  - ADF-test
  - iterative KPSS-test figure
  - Plots of the observables
  - Summary statistics
* 'StabilityScript.m' - generates the following tables/tests
  - Chow-tests
  - Cusum-tests
  - Standardized observable plots
