% run_fim_multi.m
% Multi-target OFDM FIM (no-ICI and ICI) for the reduced EVA channel (3 paths).
% Reduced EVA tap profile (Table 2 of eva_channel.tex / Section 3.1 of init):
%   tau = [0, 370, 150] ns,  rel. power = [0, -0.6, -1.4] dB
% Powers are renormalized so that sum |alpha_p|^2 = 1.

clear; clc;
addpath('../crb/EVA');

% ---- Grid / channel parameters ----
N      = 512;
M      = 128;
deltaf = 15e3;
T      = 1/deltaf;
Tcp    = 0;
Ts     = T + Tcp;
sigma2 = 1;

c        = 3e8;
fc       = 4e9;
velocity = 120 * 1e3 / 3600;        % 120 km/h
nu_max   = velocity * fc / c;

% ---- Reduced EVA channel (P = 3) ----
time_delays_ns  = [0,  370,  150];
relative_powers = [0,  -0.6, -1.4];                 % dB

P     = numel(time_delays_ns);
taus  = time_delays_ns * 1e-9;                      % [s]
P_lin = 10.^(relative_powers/10);
P_lin = P_lin / sum(P_lin);                          % unit-energy normalization
amps  = sqrt(P_lin);

rng(0);
% Per-path random phase and Doppler shift (Jakes)
phases = 2*pi*rand(P,1);
thetas = 2*pi*rand(P,1);
alphas = amps(:) .* exp(1j*phases);
nus    = nu_max * cos(thetas);

% ---- Unit-energy QPSK data ----
bits = randi([0 3], N, M);
X    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));

% ---- Compute FIMs ----
J_noICI = getFIM_OFDM_multi_noICI(X, deltaf, Ts, sigma2, alphas, taus, nus);
J_ICI   = getFIM_OFDM_multi_ICI  (X, deltaf, Ts, sigma2, alphas, taus, nus);

% ---- CRBs (symmetric diagonal preconditioning to keep inv() well-scaled) ----
CRB_noICI = invFIM_scaled(J_noICI);
CRB_ICI   = invFIM_scaled(J_ICI);

% ---- Report ----
fprintf('Reduced EVA: P = %d paths\n', P);
fprintf('%-4s %-10s %-12s %-10s %-12s\n', 'p', '|alpha|', 'tau [ns]', 'nu [Hz]', '|nu/numax|');
for p = 1:P
    fprintf('%-4d %-10.4f %-12.1f %-10.2f %-12.3f\n', ...
        p, abs(alphas(p)), taus(p)*1e9, nus(p), nus(p)/nu_max);
end
fprintf('\n');

fprintf('FIM size           : %d x %d\n', size(J_noICI,1), size(J_noICI,2));
fprintf('cond(J) no-ICI     : %.3e\n', cond(J_noICI));
fprintf('cond(J)    ICI     : %.3e\n', cond(J_ICI));
fprintf('||J_ICI - J_noICI||_F / ||J_noICI||_F = %.3e\n', ...
        norm(J_ICI - J_noICI, 'fro') / norm(J_noICI, 'fro'));
fprintf('\n');

% Per-path CRBs for tau and nu
fprintf('%-4s %-15s %-15s %-15s %-15s\n', 'p', 'CRB_tau noICI', 'CRB_tau ICI', 'CRB_nu noICI', 'CRB_nu ICI');
for p = 1:P
    i_tau = 4*(p-1)+3;
    i_nu  = 4*(p-1)+4;
    fprintf('%-4d %-15.3e %-15.3e %-15.3e %-15.3e\n', p, ...
        CRB_noICI(i_tau,i_tau), CRB_ICI(i_tau,i_tau), ...
        CRB_noICI(i_nu, i_nu ), CRB_ICI(i_nu, i_nu ));
end

% ---------- local helper ----------
function C = invFIM_scaled(J)
    d  = sqrt(diag(J));
    Js = J ./ (d * d.');
    Cs = inv(Js);
    C  = Cs ./ (d * d.');
end
