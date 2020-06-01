function [FE, Fit, Obs, Dates_Test] = FC_VAR(Y, plag, StartOoS, WinSize, Dates, h)
% Forecast a VAR model.

t      = StartOoS; 
Window = t-WinSize-1:t-1;
Y_Train = Y(Window,:); 
%Y_Test  = Y(t-1+h,:);
%Dates_Test = Dates(t:t+h-1);

[VAR, VARopt] = VARmodel(Y_Train,plag,0);

Endo = [];
for i = 0:plag-1;
    %Endo = [Endo, Y(Window(end)-i,:)];
    Endo = [Endo, VAR.ENDO(end-i,:)];
end 

FT = VAR.Fcomp;
iter = 0;
for i = h
    iter   = iter + 1; 
    Fitfoo = FT^i * Endo';
    Fitfoo = Fitfoo(1:size(Y,2))';
    Fit(:,iter) = Fitfoo;
    Y_Test = Y(t-1+i,:);
    FE(:,:,iter)  = Y_Test - Fitfoo;
    Dates_Test(iter,:) = Dates(Window(end)+i);
    Obs(:,:,iter) = Y_Test;              % Observed
end
end

