function FAVARopts = GenerateOptionsFAVAR(typeRest, Conf, nsteps, vnames)
% Auxiliary function generating options for the FAVAR, for use with the
% SVAR toolbox. 
% As most of the specifications used are similar, this functions main
% purpose is keeping the main file more tidy.
% Author: Rasmus M. Jensen

if typeRest == 1
    type = 'oir';
else
    type = 'bq';
end

if Conf == 68
    nStd = 1;
else 
    nStd = 2; 
end
FAVARopts = VARoption;
FAVARopts.nsteps = nsteps;
FAVARopts.ident  = type; % 'oir' if Cholesky, 'bq' if Blanchard-Quah
FAVARopts.impact = 0 ; % 0 for 1 st.dev. shock, 1 for unitary shock
FAVARopts.HD_pickShock   = [0]; % selected shock for HD plots
FAVARopts.method = 'bs'; % 'bs' for a standard residual bootstrapping, 'wild' for wild bootstrap (robust to heteroskedasticity)
FAVARopts.pctg = Conf; %confidence level
FAVARopts.ndraws = 1000;% number of bootstrap draws
FAVARopts.suptitle=1; % just for plotting if using the native Cesa-Bianchi
FAVARopts.vnames =    vnames;
FAVARopts.nStd = nStd;
FAVARopts.NegativeShockIndex = [];
FAVARopts.CIMethod = 1; 
end