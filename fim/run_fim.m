% run_fim.m
% Compute single-target OFDM FIM (no-ICI and ICI) using the first EVA path only.
% Parameter ordering follows ofdm_fim.tex: eta = [alphaR, alphaI, tau, nu].

clear; clc;
addpath('../crb/EVA');   % access get_ici_coeffs / get_ici_coeffs_derivative

% ---- Grid / channel parameters (mirror crb/EVA/initialize_parameters.m) ----
N      = 512;             % subcarriers
M      = 128;             % OFDM symbols
deltaf = 15e3;            % subcarrier spacing [Hz]
T      = 1/deltaf;
Tcp    = 0;
Ts     = T + Tcp;
sigma2 = 1;

c        = 3e8;
fc       = 4e9;
velocity = 120 * 1e3 / 3600;        % 120 km/h
nu_max   = velocity * fc / c;       % ~444 Hz

% ---- Use ONLY the first EVA path (tau_0 = 0, P_0 = 0 dB after normalization) ----
tau   = 0;                          % first-path delay
% Pick a representative Doppler for the first path (Jakes draw); use nu_max for a strong-ICI demo.
nu    = nu_max;
alpha = 1;                          % unit-power single path (normalized)

% ---- Unit-energy QPSK data symbols ----
rng(0);
bits  = randi([0 3], N, M);
X     = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));

% ---- Compute FIMs ----
J_noICI = getFIM_OFDM_noICI(X, deltaf, Ts, sigma2, alpha, tau, nu);
J_ICI   = getFIM_OFDM_ICI (X, deltaf, Ts, sigma2, alpha, tau, nu);

names = {'alphaR','alphaI','tau','nu'};

fprintf('FIM (no-ICI):\n');     disp(array2table(J_noICI,'VariableNames',names,'RowNames',names));
fprintf('FIM (ICI):\n');        disp(array2table(J_ICI  ,'VariableNames',names,'RowNames',names));

% CRBs for tau and nu (diagonal of inverse FIM)
CRB_noICI = inv(J_noICI);
CRB_ICI   = inv(J_ICI);

fprintf('CRB(tau)  no-ICI = %.3e   ICI = %.3e   [s^2]\n', CRB_noICI(3,3), CRB_ICI(3,3));
fprintf('CRB(nu)   no-ICI = %.3e   ICI = %.3e   [Hz^2]\n', CRB_noICI(4,4), CRB_ICI(4,4));
