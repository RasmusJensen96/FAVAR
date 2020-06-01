function [Ytph_Filtered phi FE OBSAR1] = AR_FC(Y, t, WinSize, h)
% Forecast an Autoregressive model of order 1 to provided data vector by OLS
% 
% Rasmus M. Jensen 2020

Window = t-WinSize-1:t-1;
Ytest = Y(Window,:); 
Yt = Ytest(4:end); %Ytm1 = Ytest(1:end-1);
YtL = [Ytest(3:end-1), Ytest(2:end-2), Ytest(1:end-3)];
%Yt   = Y(h:t);
%Ytm1 = Y(1:t-1);


Ytph = Y(t+h-1);
% phi = inv(Ytm1'*Ytm1) * Ytm1'*Yt;
%phi = nwest(Yt, Ytm1,1);
phi = nwest(Yt, YtL,1);
phi = phi.beta;

% One period ahead forecast:
%Ytph_Filtered = Yt(end) * (phi')^h;
foodata = Yt(end-2:end);
for i = 1:h
    foo = foodata' * phi;
    foodata = [foo; foodata(2:end)];
end
Ytph_Filtered =   foo;
%Ytph_Filtered = Yt(end-3:end)' * phi^h;
OBSAR1 = Ytph;
FE = Ytph-Ytph_Filtered;
end