function [LLH, Sigma2] = GARCHFilter(vY,dOmega, dAlpha, dBeta)
% Soft enforced constraint
foo = dAlpha + dBeta;
if foo > 1 || foo<0 || dAlpha>1 || dBeta>1
      LLH = -1e9; 
      Sigma2 = [];
      return
end

  iT = size(vY,1);
  Sigma2 = zeros(iT,1);
  Sigma2(1) = dOmega/(1-dAlpha-dBeta); %var(vY)
  %Sigma2(1) = var(vY(1:round(iT * 0.1)));
  %Sigma2(1) = 1; % for our VAR variance is initialized at 1
  
  % compute the log--likelihood of the first obs
  LLH = log(normpdf(vY(1), 0, sqrt(Sigma2(1))));
  % loop over iT
  for t = 2:iT 
    %  update the volatility
    Sigma2(t) = dOmega + dAlpha * vY(t-1)^2 + dBeta * Sigma2(t - 1);
    % add the log-likelihood contribution at time t
    LLH = LLH + log(normpdf(vY(t), 0, sqrt(Sigma2(t))));
  end
end