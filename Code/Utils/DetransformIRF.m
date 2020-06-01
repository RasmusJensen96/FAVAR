function [IRF, Upper, Lower] = DetransformIRF(uaIRF, uaUIRFR, uaLIRF, tcodes)
% =========================================================================
% DESCRIPTION: 
% This function untransforms transformed IRFs based on each series' transformation
% code.
%
% -------------------------------------------------------------------------
% INPUT:
%           IRFs and conf inv.     = raw data 
%           tcode       = transformation codes for each series
%
% OUTPUT: 
%           yt          = transformed data
%
% -------------------------------------------------------------------------
% SUBFUNCTION:
%           untransxf:    transforms a single series as specified by a 
%                       given transfromation code
%
% =========================================================================
% APPLY TRANSFORMATION:
% Initialize output variable
IRF        = [];
Upper      = [];
Lower      = [];

% Number of series kept
N = size(uaIRF,1);                         

% Perform transformation using subfunction transxf (see below for details)
for i = 1:N
    [Med1, Upper1, Lower1] = untransxf(uaIRF(i,:), uaUIRFR(i,:), uaLIRF(i,:),tcodes(i));
    IRF    = [IRF, Med1']; Upper = [Upper, Upper1']; Lower = [Lower, Lower1'];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION

function [Med1, Upper1, Lower1] =untransxf(uaIRF, uaUIRFR, uaLIRF,tcodes)
% =========================================================================
% DESCRIPTION:
% This function transforms a single series (in a column vector)as specified
% by a given transfromation code.
%
% -------------------------------------------------------------------------
% INPUT:
%           x       = series (in a column vector) to be transformed
%           tcode   = transformation code (1-7)
%
% OUTPUT:   
%           Med upper lower = transformed series (as a column vectors)
%
% =========================================================================
% TRANSFORMATION: 
% Determine case 1-7 by transformation code
n = size(uaUIRFR, 2);
switch(tcodes)
  case 1 % Level (i.e. no transformation): x(t)
    Med1=uaIRF;
    Lower1 = uaLIRF;
    Upper1 = uaUIRFR;

  case 2 % First difference: x(t)-x(t-1)
    Med1(1:n)   = cumsum(uaIRF);
    Lower1(1:n) = cumsum(uaLIRF);
    Upper1(1:n) = cumsum(uaUIRFR);
  
  case 3 % Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
    Med1(1:n)   = cumsum(cumsum(uaIRF));
    Lower1(1:n) = cumsum(cumsum(uaLIRF));
    Upper1(1:n) = cumsum(cumsum(uaUIRFR));

  case 4 % Natural log: ln(x)
    Med1(1:n)   = exp(uaIRF);
    Lower1(1:n) = exp(uaLIRF);
    Upper1(1:n) = exp(uaUIRFR);
  
  case 5 % First difference of natural log: ln(x)-ln(x-1)
    Med1(1:n)   = exp(cumsum(uaIRF))-1;
    Lower1(1:n) = exp(cumsum(uaLIRF))-1;
    Upper1(1:n) = exp(cumsum(uaUIRFR))-1;
  
  case 6 % Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
    %Med1(1:n)    = exp(cumsum(cumsum(uaIRF)))-1;
    %Lower1(1:n)  = exp(cumsum(cumsum(uaLIRF)))-1;
    %Upper1(1:n)  = exp(cumsum(cumsum(uaUIRFR)))-1;
    Med1(1:n)    = [exp(uaIRF(1))-1, exp(cumsum(uaIRF(2:n)) + cumsum(uaIRF(1:n-1)))-1];% exp(cumsum(uaIRF(2:n+1)) + cumsum(uaIRF(1:n)));
    Lower1(1:n)  = [exp(uaLIRF(1))-1, exp(cumsum(uaLIRF(2:n)) + cumsum(uaLIRF(1:n-1)))-1];
    Upper1(1:n)  = [exp(uaUIRFR(1))-1, exp(cumsum(uaUIRFR(2:n)) + cumsum(uaUIRFR(1:n-1)))-1];
  case 7
    Med1(1:n)   = exp(cumsum(uaIRF))-1;
    Lower1(1:n) = exp(cumsum(uaLIRF))-1;
    Upper1(1:n) = exp(cumsum(uaUIRFR))-1;
end
end