% [nfact, v_nfact, cr] = NumFacsDFM_FAVAR(panel, q_max, band, nbck, stp, c_max, penalty, cf, m, h, plot_opt)
% 
% log criterion to determine the number of dynamic factors according to 
% Hallin and Liska (2007) "Determining the Number of Factors in the General 
% Dynamic Factor Model", Journal of the American Statistical Association, 
% 102, 603-617 
%
% Original function: Hallin and Liska (2007), name (numfactors_nonstd_log_band.m)
% Modified:          Rasmus M. Jensen (2020)
%
% INPUT:    panel           :   T x n data matrix 
%                               data should be covariance stationary 
%           q_max           :   upper bound on the number of factors
%           band            :   restrict to a frequency band 
%                               [-band(2)*pi : -band(1)*pi , band(1)*pi : band(2)*pi]
%                               (default value: 1 1)
%           nbck, stp       :   T x n_j subpanels are used where
%                               n_j = n - nbck : stp: n 
%                               (default value: nbck = floor(n/4), stp = 1)
%           c_max           :   c = [0:cmax] (default value: 3)
%           penalty         :   p1 = ((m/T)^0.5 + m^(-2) + n^(-1))*log(min([(T/m)^0.5;  m^2; n]))  
%                               p2 = (min([(T/m)^0.5;  m^2; n])).^(-1/2)  
%                               p3 = (min([(T/m)^0.5;  m^2; n])).^(-1)*log(min([(T/m)^0.5;  m^2; n]))
%                               (default value: 'p1')
%           cf              :   1/cf is granularity of c 
%                               (default value: 1000)
%           m               :   covariogram truncation 
%                               (default value: floor(sqrt(T)))
%           h               :   number of points in which the spectral 
%                               density is computed (default value: m)
%           plot_opt        :   option to draw the plot 
%                               (yes == 1, no == 0)(default value: 1)
%
% OUTPUT:   nfact           :   number of dynamic factors as function of c
%                               computed for n_j = n
%           v_nfact         :   variance in the number of dynamic factors
%                               as function of c and computed as the 
%                               n_j varies  
%           cr              :   values of c (needed for the plot)

function [nfact, v_nfact, cr] = NumFacsDFM_FAVAR(panel, q_max, band, nbck, stp, c_max, penalty, cf, m, h, plot_opt)

%% Preliminary settings
[T,n] = size(panel);

if nargin < 2 
    disp('ERROR MESSAGE: Too few input arguments'); 
    return 
end

if nargin==2
%     band=[1 1];
    band=[1 1];
    nbck = floor(n/4);
    stp = 1;
    c_max = 5;
    penalty = 'p2';
    cf = 1000;
    m = floor(sqrt(T));
    h = m;
    plot_opt = 1;
end

if nargin == 3
    nbck = floor(n/4);
    stp = 1;
    c_max = 3;
    penalty = 'p1';
    cf = 1000;
    m = floor(sqrt(T));
    h = m;
    plot_opt = 1;
end

if nargin == 4
    stp = 1;
    c_max = 3;
    penalty = 'p1';
    cf = 1000;
    m = floor(sqrt(T));
    h = m;
    plot_opt = 1;
end

if nargin == 5
    c_max = 3;
    penalty = 'p1';
    cf = 1000;
    m = floor(sqrt(T));
    h = m;
    plot_opt = 1;
end

if nargin == 6   
    penalty = 'p1';
    cf = 1000;
    m = floor(sqrt(T));
    h = m;
    plot_opt = 1;
end

if strcmp(penalty, 'p1') == 0 && strcmp(penalty, 'p2') == 0 && strcmp(penalty, 'p3') == 0
    disp('ERROR MESSAGE : Penalty function can only take 3 values: p1, p2 and p3');
    return
end

if nargin == 7
    cf = 1000;
    m = floor(sqrt(T));
    h = m;
    plot_opt = 1;
end

if nargin == 8
    m = floor(sqrt(T));
    h = m;
    plot_opt = 1;
end

if nargin == 9
    h = m;
    plot_opt = 1;
end

if nargin == 10
    plot_opt = 1;
end

%% Mean-standardize data
m_X = mean(panel);
%s_X = std(panel);
X = (panel - ones(T,1)*m_X); %./s_X;

%% Compute the number of dynamic factors
s=0;
for N = n-nbck:stp:n
    disp(sprintf('subsample size %d',N));
    s = s+1;
    [a rv] = sort(rand(n,1));                                                 % select randomly N series
    subpanel = X(1:T,rv(1:N));

    m_subpanel = mean(subpanel);
%     s_subpanel = std(subpanel);
    subpanel = (subpanel - ones(T,1)*m_subpanel);                           % standardize the subpanel

    [P_X, D_X, Sigma_X] = spectral_band(subpanel, N, h, m, band);           % in this case we use spectral with q = N
%     E = [D_X(:,h+1)  D_X(:,h+2:2*h+1)*2]*ones(h+1,1)/(2*h+1);             % all the n dynamic eigenvalues
    E = [D_X(:,h+1:2*h)*2]*ones(h,1)/(2*h);                                 % all the n dynamic eigenvalues
    IC1 = flipud(cumsum(flipud(E)));                                        % compute information criterion
    IC1 = IC1(1:q_max+1,:);
    
    if strcmp(penalty, 'p1') == 1                                   
        p = ((m/T)^0.5 + m^(-2) + N^(-1))*log(min([(T/m)^0.5;  m^2; N]))*ones(q_max+1,1);  
    elseif strcmp(penalty, 'p2') == 1
        p = (min([(T/m)^0.5;  m^2; N])).^(-1/2)*ones(q_max+1,1);  
    elseif strcmp(penalty, 'p3') == 1    
        p = (min([(T/m)^0.5;  m^2; N])).^(-1)*log(min([(T/m)^0.5;  m^2; N]))*ones(q_max+1,1);  
    end

    for c = 1:floor(c_max*cf)
        cc = c/cf;
        IC = log(IC1./N) + (0:q_max)'.*p*cc;
        rr = find((IC == ones(q_max+1,1)*min(IC))==1);              % compute minimum of IC
        o(s,c) = rr-1;
    end
end

cr = (1:floor(c_max*cf))'/cf;
nfact = o(end,:);                                                       % number of factors when N = n
v_nfact = std(o);

%% Plot if needed
if plot_opt == 1
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    figure
    yyaxis left
    plot(cr,nfact,'k-')
    ylabel('Number of dynamic factors','interpreter','latex')
    yAX = get(gca,'YAxis');yAX(1).FontSize = 18;yAX(1).TickLabelInterpreter = "latex";
    yAX(2).FontSize = 18;yAX(2).TickLabelInterpreter = "latex";
    xAX = get(gca,'XAxis');xAX.FontSize = 18;xAX.TickLabelInterpreter = "latex";
    hold on
    yyaxis right
    plot(cr,v_nfact,'k:')
    ylabel('Standard-deviation of $\hat{K}$ as $n \rightarrow N$','interpreter','latex')
    hold all
    xlabel('Penalty Weight, $c$','interpreter','latex')
    %yAX = get(gca,'YAxis');yAX.FontSize = 18;yAX.TickLabelInterpreter = "latex";
    %xAX = get(gca,'XAxis');xAX.FontSize = 18;xAX.TickLabelInterpreter = "latex";
    legend('$S_c$','$K$','interpreter','latex','FontSize',18)
    set(gcf, 'PaperPosition', [0 0 40 30]);
    set(gcf, 'PaperSize', [40 30]);
    %title('estimated number of factors - ','interpreter','latex')
end
end

function w = triang(n_est)
%TRIANG Triangular window.
%   W = TRIANG(N) returns the N-point triangular window.
%
%   See also BARTHANNWIN, BARTLETT, BLACKMANHARRIS, BOHMANWIN, 
%            FLATTOPWIN, NUTTALLWIN, PARZENWIN, RECTWIN, WINDOW.

%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.13.4.1 $  $Date: 2007/12/14 15:06:30 $

error(nargchk(1,1,nargin,'struct'));
[n,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

if rem(n,2)
    % It's an odd length sequence
    w = 2*(1:(n+1)/2)/(n+1);
    w = [w w((n-1)/2:-1:1)]';
else
    % It's even
    w = (2*(1:(n+1)/2)-1)/n;
    w = [w w(n/2:-1:1)]';
end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [P_chi, D_chi, Sigma_chi] = spectral(X, q, m, h, band)
%
% Computes spectral decomposition for the data matrix in input.
%
% INPUT:    X               :   T x n data matrix
%                               data should be covariance stationary and
%                               mean standardized
%           q               :   dimension of the common space
%                               (i.e. number of dynamic factors)
%           m               :   covariogram truncation
%                               (default value: floor(sqrt(T)))
%           h               :   number of points in which the spectral
%                               density is computed (default value: m)
%           band            :   restrict to a frequency band 
%                               [-band(2)*pi : -band(1)*pi , band(1)*pi : band(2)*pi]
%           
% 
%
% OUTPUT:   P_chi           :   n x q x 2*h+1 matrix of dynamic eigenvectors
%                               associated with the q largest dynamic eigenvalues
%                               for different frequency levels
%           D_chi           :   q x 2*h+1 matrix of dynamic eigenvalues
%                               for different frequency levels
%           Sigma_chi       :   (n x n x 2*h+1) spectral density matrix
%                               of common components with the 2*h+1 density
%                               matrices for different frequency levels
function [P_chi, D_chi, Sigma_chi] = spectral_band(X, q, m, h, band)

%% Preliminary settings
[T,n] = size(X);

if nargin < 2
    disp('ERROR MESSAGE: Too few input arguments');
    return
end

if nargin == 2
    m = floor(sqrt(T));
    h = m;
    band=[1 1];
end

if nargin == 3
    h = m;
%     band=[1 1];
    band=[];
end

%% Compute M covariances
M = 2*m+1;
B = triang(M);                                                              % Triangular window (similar Bartlett)
Gamma_k = zeros(n,n,M);
for k = 1:m+1,
    Gamma_k(:,:,m+k) = B(m+k)*(X(k:T,:))'*(X(1:T+1-k,:))/(T-k);
    Gamma_k(:,:,m-k+2) = Gamma_k(:,:,m+k)';
end

%% Compute the spectral density matrix in H points
% if all(band==[1 1])
if isempty(band)
    H = 2*h+1;
    Factor = exp(-sqrt(-1)*(-m:m)'*(-2*pi*h/H:2*pi/H:2*pi*h/H));
else
    H = 2*h;
    freq=linspace(band(1)*2*pi*h/H,band(2)*2*pi*h/H,h);
    freq=[-freq freq(end:-1:1)];
    Factor = exp(-sqrt(-1)*(-m:m)'*freq);
end


Sigma_X = zeros(n,n,H);
for j = 1:n
    Sigma_X(j,:,:) = squeeze(Gamma_k(j,:,:))*Factor;
end

%% Create output elements
P_chi = zeros(n,q,H); 
D_chi = zeros(q,H);
Sigma_chi = zeros(n,n,H);

%% Compute eigenvalues and eigenvectors
%% case with q < n-1 we can use eigs, fast method
if q < n-1 
    opt.disp = 0;
    if all(band==[1 1])
        [P, D] = eigs(squeeze(Sigma_X(:,:,h+1)),q,'LM',opt);                    % frequency zero
        D_chi(:,h+1) = diag(D);
        P_chi(:,:,h+1) = P;
        Sigma_chi(:,:,h+1) = P*D*P';
    end
    
    for j = 1:h
        [P, D] = eigs(squeeze(Sigma_X(:,:,j)),q,'LM',opt);                  % other frequencies
        D_chi(:,j) = diag(D);
        D_chi(:,H+1-j) = diag(D);
        P_chi(:,:,j) = P;
        P_chi(:,:,H+1-j) = conj(P);

        Sigma_chi(:,:,j) = P*D*P';
        Sigma_chi(:,:,H+1-j) = conj(P*D*P');
    end
end

%% case with q >= n-1, we must use eig, slow method
if q >= n-1 
    [P, D] = eig(squeeze(Sigma_X(:,:,h+1)));                                % frequency zero
    [D,IX] = sort((diag(D)));                                               % sort eigenvalues and eigenvectors
    D = flipud(D);
    IX = flipud(IX);
    P = P(:,IX);
    D = diag(D);
    D_chi(:,h+1) = real(diag(D));
    P_chi(:,:,h+1) = P;
    Sigma_chi(:,:,h+1) = P*D*P';

    for j = 1:h
        [P, D] = eig(squeeze(Sigma_X(:,:,j)));                              % other frequencies
        [D,IX] = sort((diag(D)));                                           % sort eigenvalues and eigenvectors
        D = flipud(D);
        IX = flipud(IX);
        P = P(:,IX);
        D = diag(D);
        D_chi(:,j) = real(diag(D));
        D_chi(:,H+1-j) = real(diag(D));
        P_chi(:,:,j) = P;
        Sigma_chi(:,:,j) = P*D*P';
        P_chi(:,:,H+1-j) = conj(P);
        Sigma_chi(:,:,H+1-j) = conj(P*D*P');
    end 
end
end