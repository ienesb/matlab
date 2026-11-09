% run_fim_otfs.m
% Single-target OTFS FIM using the first EVA path only.
% Parameter ordering follows otfs_fim.tex: eta = [alphaR, alphaI, tau, nu].

clear; clc;
addpath('../crb/EVA');

% ---- Grid / channel parameters (mirror crb/EVA/initialize_parameters.m) ----
N      = 512;
M      = 128;
deltaf = 15e3;
T      = 1/deltaf;
Tcp    = 0;
Ts     = T + Tcp;
sigma2 = 1;

c        = 3e8;
fc       = 4e9;
velocity = 120 * 1e3 / 3600;
nu_max   = velocity * fc / c;

% ---- First EVA path (tau_0 = 0, unit power for a single-path channel) ----
tau   = 0;
nu    = nu_max;
alpha = 1;

% ---- Unit-energy QPSK data in the delay-Doppler grid ----
rng(0);
bits = randi([0 3], N, M);
x    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));

% ---- OTFS FIM ----
J_OTFS = getFIM_OTFS(x, deltaf, Ts, sigma2, alpha, tau, nu);

% ---- Cross-check: OFDM no-ICI FIM with the same QPSK symbols (placed in TF) ----
% otfs_fim.tex (paragraph after eq. for E{Jnn}) proves the two expected FIMs
% are identical entry-by-entry under i.i.d. unit-energy QPSK. With a single
% finite realization the match should be very close, not exact.
J_OFDM_noICI = getFIM_OFDM_noICI(x, deltaf, Ts, sigma2, alpha, tau, nu);

names = {'alphaR','alphaI','tau','nu'};
fprintf('OTFS FIM:\n');           disp(array2table(J_OTFS,      'VariableNames',names,'RowNames',names));
fprintf('OFDM no-ICI FIM (same QPSK realization):\n');
disp(array2table(J_OFDM_noICI, 'VariableNames',names,'RowNames',names));

fprintf('||J_OTFS - J_OFDMnoICI||_F / ||J_OFDMnoICI||_F = %.3e\n', ...
        norm(J_OTFS - J_OFDM_noICI, 'fro') / norm(J_OFDM_noICI, 'fro'));

% CRBs for tau and nu (use scaled inverse to dodge the ill-scaled diagonals)
CRB_OTFS      = invFIM_scaled(J_OTFS);
CRB_OFDMnoICI = invFIM_scaled(J_OFDM_noICI);

fprintf('\n');
fprintf('CRB(tau): OTFS = %.3e   OFDM no-ICI = %.3e   [s^2]\n', CRB_OTFS(3,3), CRB_OFDMnoICI(3,3));
fprintf('CRB(nu) : OTFS = %.3e   OFDM no-ICI = %.3e   [Hz^2]\n', CRB_OTFS(4,4), CRB_OFDMnoICI(4,4));

% ---------- local helper ----------
function C = invFIM_scaled(J)
    d  = sqrt(diag(J));
    Js = J ./ (d * d.');
    Cs = inv(Js);
    C  = Cs ./ (d * d.');
end
