%% Call TFAVAR
tic;
addpath(genpath('Data'));
addpath(genpath('Utils'));
clear all; clc;
set(0,'defaultAxesFontSize',14,'defaultTextInterpreter','latex');
rng(260196,'twister');
MLGARCH          = 1;
Cholesky         = 0;
signRest         = 0;  % Aux switch indicating whether or not sign-restrictions should be applied
faclag           = 3;
FFR1d            = 1;  % Switch to determine whether FFR should be assumed I(1) or levels I(0)
type             = 2;  % Switch 1 = PCA extraction of factors, 2 = Kalman Smoothing of factors. NB. Kalman is computationally very heavy.
nsteps           = 24; % Number of steps IRF
ConfLevel        = 68; % Confidence level percentage for bootstrap
if Cholesky == 1
    Order = [1,2,3];
else
    Order = [2,3,1];
end
%%
PrepData_TFAVAR;
if FFR1d == 1
%    ydataPre_st(:,3) = [NaN; ydataPre_st(2:end,3)-ydataPre_st(1:end-1,3)];
%    ydataPre_st = ydataPre_st(3:end,:);
%    ydataPost_st(:,3) = [NaN; ydataPost_st(2:end,3)-ydataPost_st(1:end-1,3)];
%    ydataPost_st = ydataPost_st(3:end,:);
    firstdstr = '1d';
%    plag = 3;
else
 %   ydataPre_st = ydataPre_st(3:end,:);
 %   ydataPost_st = ydataPost_st(3:end,:);
    firstdstr = 'Level';
 %   plag = 3;
end
Xpre = x2Pre_st;
plag = 3;
K = 3;
det = 0;
Xpost = x2Post_st;
%%
[FAVARPre, FAVARPreopt]   = EstimateReducedFAVAR(Xpre, ydataPre_st, faclag, plag, K, det, slowcode, 1,namesXY,MLGARCH,Cholesky);
[FAVARPost, FAVARPostopt] = EstimateReducedFAVAR(Xpost, ydataPost_st, faclag, plag, K, det, slowcode, type,namesXY,MLGARCH,Cholesky);
[FAVARPreopt.standardized, FAVARPostopt.standardized] = deal(1);
%%
vnames = ["Factor 1", "Factor 2", "Factor 3", namesY'];
if Cholesky == 1
FAVARCholopt = GenerateOptionsFAVAR(1, ConfLevel, nsteps,vnames);
% Compute IRF
[IRFCholPre, FAVARCholPre] = VARir(FAVARPre,FAVARCholopt);
[IRFCholPost, FAVARCholPost] = VARir(FAVARPost,FAVARCholopt);
% Compute error bands
[IRFINF_Chol,IRFSUP_Chol,IRFMED_Chol, INFxChol, SUPxChol, MEDxChol] = FAVARirband(FAVARCholPre,FAVARCholopt);
save("IRFMedCholPre", 'IRFMED_Chol')
[IRFINF_CholPost,IRFSUP_CholPost,IRFMED_CholPost, INFxCholPost, SUPxCholPost, MEDxCholPost] = FAVARirband(FAVARCholPost,FAVARCholopt);
% Plot 
nlagstr = num2str(FAVARCholPre.nlag);
name  = ['chol',nlagstr,'LagsThresPreSample',firstdstr];
FAVAR_Resp(FAVARCholPre, FAVARCholopt, IRFMED_Chol, IRFINF_Chol, IRFSUP_Chol,name, FFR1d,false);
name  = ['chol',nlagstr,'LagsThresPostSample',firstdstr];
FAVAR_Resp(FAVARCholPost, FAVARCholopt, IRFMED_CholPost, IRFINF_CholPost, IRFSUP_CholPost,name, FFR1d,false);
else
FAVARbqopt = GenerateOptionsFAVAR(2, ConfLevel, nsteps,vnames);
FAVARbqopt.standardized = 0;
% Compute IRF
[IRFbqPre, FAVARPre] = VARir(FAVARPre,FAVARbqopt);
[IRFbqPost, FAVARPpost] = VARir(FAVARPost,FAVARbqopt);
% Compute error bands
[IRFINFPre_bq,IRFSUPPre_bq,IRFMEDPre_bq, INFxPre, SUPxPre, MEDxPre] = FAVARirband(FAVARPre,FAVARbqopt);
[IRFINFPost_bq,IRFSUPPost_bq,IRFMEDPost_bq, INFxPost, SUPxPost, MEDxPost] = FAVARirband(FAVARPost,FAVARbqopt);
save("IRFMedbqPre", 'IRFMEDPre_bq')
% Plot 
nlagstr = num2str(FAVARPre.nlag);
name  = ['bq',nlagstr,'LagsThresPreSample',firstdstr];
FAVAR_Resp(FAVARPre,  FAVARbqopt, IRFMEDPre_bq, IRFINFPre_bq, IRFSUPPre_bq,name, FFR1d,false);
%%
name  = ['bq',nlagstr,'LagsThresPostSample',firstdstr];
FAVAR_Resp(FAVARPost, FAVARbqopt, IRFMEDPost_bq, IRFINFPost_bq, IRFSUPPost_bq,name, FFR1d,false);
end
%%
TotFit = [FAVARPre.fit; FAVARPost.fit];
TotRes = [FAVARPre.residuals; FAVARPost.residuals];
sqrt(mean(TotRes.^2))
%%
figure
stackedplot(TotRes); title('Residuals');
figure
stackedplot(TotFit); title('Filtered')