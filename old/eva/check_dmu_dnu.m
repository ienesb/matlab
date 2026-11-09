% check_dmu_dnu.m
% Finite-difference test of the Doppler derivative dmu/dnu used by the ICI
% FIM functions, against the updated ofdm_fim.pdf (eq. 8 with I' = (1/j2pi) dI/dnu).
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

deltaf = 15e3;  T = 1/deltaf;  Ts = T;
N = 8; M = 5;
alpha = 0.8*exp(1j*0.7);  tau = 370e-9;  nu = 300;

rng(1);
b = randi([0 3],N,M);
X = (1/sqrt(2))*((2*bitand(b,1)-1)+1j*(2*floor(b/2)-1));

mu = @(v) mu_ici(X, deltaf, Ts, alpha, tau, v, N, M);

% Central finite difference (truth).
h = 1e-3;
dmu_fd = (mu(nu+h) - mu(nu-h)) / (2*h);

% Analytic, two candidate conventions, both using kernel + dI/dnu.
ns = (0:N-1).'; ms = (0:M-1);
[NP,NN] = ndgrid(0:N-1,0:N-1); diff_n = NP-NN;
I_mat  = get_ici_coeffs(diff_n, nu, T);
dI_mat = get_ici_coeffs_derivative(diff_n, nu, T);   % = dI/dnu (full)
A = X .* exp(-1j*2*pi*deltaf*ns*tau) .* exp(1j*2*pi*Ts*ms*nu);
Phi   = I_mat.'  * A;
dIcon = dI_mat.' * A;

% (A) CODE convention: g_nu = j2pi*alpha*(mTs*Phi + dI/dnu term)
gnu_code = 1j*2*pi*alpha * ((ms*Ts).*Phi + dIcon);
% (B) TRUE / PDF: dmu/dnu = alpha*(j2pi*mTs*Phi + dI/dnu term)
gnu_true = alpha * (1j*2*pi*(ms*Ts).*Phi + dIcon);

fprintf('||code - FD|| / ||FD|| = %.3e\n', norm(gnu_code(:)-dmu_fd(:))/norm(dmu_fd(:)));
fprintf('||true - FD|| / ||FD|| = %.3e\n', norm(gnu_true(:)-dmu_fd(:))/norm(dmu_fd(:)));

function out = mu_ici(X, deltaf, Ts, alpha, tau, nu, N, M)
    T = 1/deltaf;
    ns = (0:N-1).'; ms = (0:M-1);
    [NP,NN] = ndgrid(0:N-1,0:N-1); diff_n = NP-NN;
    I_mat = get_ici_coeffs(diff_n, nu, T);
    A = X .* exp(-1j*2*pi*deltaf*ns*tau) .* exp(1j*2*pi*Ts*ms*nu);
    out = alpha * (I_mat.' * A);   % mu[n,m]
end
