function [DCC] = Estimate_DCC(FAVAR)
%{ 
Function fitting the DCC-model to a given matrix (2 dimensional array)
    of series. (This implementation uses GARCH(1,1)-marginals)
-------------------------------------Inputs--------------------------------
            FAVAR: Structure from VARmodel
------------------------------------Outputs--------------------------------
            dLLK:  total likelihood at numerical optimum (L_v+L_c) L_v = Vol. (L_c) component 1 (2) in equation 13 in the PDF
            mCoef: Parameters of the univariate GARCH-marginals (Nx3)-dimensional (omega_j, kappa_j, lambda_j)
            vPar:  Parameters of the DCC-model (1x2)-dim (alpha, beta)
            mEta:  Standardised residuals
            BIC:   Average BIC evaluated at the numerical maximum likelihood
  
 Rasmus Jensen 2020
%}

% list where marginal models are stored
 % Init. alloacation of Eta and sigma
[iT, iN] = size(FAVAR.residuals);
mEta            = NaN(iT, iN);
mSigma          = NaN(iT, iN);

for j = 1:iN
    aux = ['eq' num2str(j)];
    mEta(:,j)   = eval(['FAVAR.' aux '.standardizedresidual;'] );  % eta
    mSigma(:,j) = eval(['FAVAR.' aux '.filteredvolatility;']); % Volatility
    LLH(j)      = eval(['FAVAR.' aux '.LLH;']); % LLH
end
LLH = -LLH;
dLLK_V = sum(LLH);
%% maximization of the DCC likelihood
  %initial parameters
  vPar = [0.05, 0.90];

  %unconditional correlation
  mQ = corr(mEta);

  %maximize the DCC likelihood
  [vPar, dLLK_C] = fminunc(@(vPar)ObjectiveDCC(vPar),vPar);
  
function [LLH] = ObjectiveDCC(vPar)
[LLH,~] = DCCFilter(mEta, vPar(1), vPar(2), mQ);
LLH = -LLH;
end

dLLK_C = -dLLK_C;

  %Filter the dynamic correlation using the estimated parameters
  [~,Rcorr] = DCCFilter(mEta, vPar(1), vPar(2), mQ);

  %compute the total likelihood
  dLLK = dLLK_V + dLLK_C;
  
for t = 1:iT
    VCV(:,:,t) = diag(mSigma(t,:)) * squeeze(Rcorr(:,:,t)) * diag(mSigma(t,:));
end
  
% Returns a list of series and values of interest.
DCC = struct();
DCC.LLH = dLLK;
DCC.Corr = Rcorr;
DCC.Sigma = mSigma;
DCC.VCV = VCV;
DCC.vPar = vPar;
DCC.eta = mEta;
end

function [LLH,Rcorr] = DCCFilter(mEta, dA, dB, mQ)
  [iT, iN] = size(mEta);
  % initialize the array for the correlations
  aCor = zeros(iN, iN, iT);
  aQ   = zeros(iN, iN, iT);
  %% initialization at the unconditional cor
  aCor(:,:,1) = mQ;
  aQ(:,:,1) = mQ;
  % Compute the first likelihood contribution
  dLLK = mEta(1,:) * inv(squeeze(aCor(:,:,1))) * mEta(1,:)' - mEta(1,:) * mEta(1,:)' + log(det(squeeze(aCor(:,:,1))));
  %main loop
  for t = 2:iT
    %update the Q matrix (Updating, GARCH(1,1), style equation for Q_t)
    aQ(:,:,t) = mQ * (1 - dA - dB) + dA * mEta(t - 1,:)' * mEta(t - 1,:) + dB * aQ(:,:,t - 1);
    %% Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2}
    foo1 = diag(diag(sqrt(1./aQ(:,:,t))));
    aCor(:,:,t) = foo1 * aQ(:,:,t) * foo1;
    %aCor(:,:,t) = diag(sqrt(1/diag(aQ(:,:,t)))) * aQ(:,:,t) * diag(sqrt(1/diag(aQ(:,:,t))));
    %augment the likelihood value
    dLLK = dLLK + mEta(t,:) * inv(squeeze(aCor(:,:,t))) * mEta(t,:)' - mEta(t,:) * squeeze(mEta(t,:,:))' + log(det(squeeze(aCor(:,:,t))));
  end
  LLH = -0.5 * dLLK;
  Rcorr = aCor;
end