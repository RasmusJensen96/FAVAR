% function [P_chi, D_chi, Sigma_chi] = spectral(X, q, m, h)
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
%
% OUTPUT:   P_chi           :   n x q x 2*h+1 matrix of dynamic eigenvectors
%                               associated with the q largest dynamic eigenvalues
%                               for different frequency levels
%           D_chi           :   q x 2*h+1 matrix of dynamic eigenvalues
%                               for different frequency levels
%           Sigma_chi       :   (n x n x 2*h+1) spectral density matrix
%                               of common components with the 2*h+1 density
%                               matrices for different frequency levels

function [P_chi, D_chi, Sigma_chi] = spectral(X, q, m, h)

%% Preliminary settings
[T,n] = size(X);

if nargin < 2
    disp('ERROR MESSAGE: Too few input arguments');
    return
end

if nargin == 2
    m = floor(sqrt(T));
    h = m;
end

if nargin == 3
    h = m;
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
H = 2*h+1;
Factor = exp(-sqrt(-1)*(-m:m)'*(-2*pi*h/H:2*pi/H:2*pi*h/H));
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
    [P, D] = eigs(squeeze(Sigma_X(:,:,h+1)),q,'LM',opt);                    % frequency zero
    D_chi(:,h+1) = diag(D);
    P_chi(:,:,h+1) = P;
    Sigma_chi(:,:,h+1) = P*D*P';

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
