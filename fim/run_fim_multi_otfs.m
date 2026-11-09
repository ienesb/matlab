% run_fim_multi_otfs.m
% Multi-target OTFS FIM for the reduced EVA channel (3 paths).
% Mirror of run_fim_multi.m for the OTFS modulation.

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
velocity = 120 * 1e3 / 3600;
nu_max   = velocity * fc / c;

% ---- Reduced EVA channel (P = 3) ----
time_delays_ns  = [0,  370,  150];
relative_powers = [0,  -0.6, -1.4];

P     = numel(time_delays_ns);
taus  = time_delays_ns * 1e-9;
P_lin = 10.^(relative_powers/10);
P_lin = P_lin / sum(P_lin);
amps  = sqrt(P_lin);

rng(0);
phases = 2*pi*rand(P,1);
thetas = 2*pi*rand(P,1);
alphas = amps(:) .* exp(1j*phases);
nus    = nu_max * cos(thetas);

% ---- Unit-energy QPSK data in the delay-Doppler grid ----
bits = randi([0 3], N, M);
x    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));

% ---- Compute FIMs (OTFS + reference OFDM no-ICI for the equivalence check) ----
J_OTFS = getFIM_OTFS_multi      (x, deltaf, Ts, sigma2, alphas, taus, nus);
J_OFDM = getFIM_OFDM_multi_noICI(x, deltaf, Ts, sigma2, alphas, taus, nus);

% ---- CRBs (symmetric diagonal preconditioning) ----
CRB_OTFS = invFIM_scaled(J_OTFS);
CRB_OFDM = invFIM_scaled(J_OFDM);

% ---- Report ----
fprintf('Reduced EVA: P = %d paths\n', P);
fprintf('%-4s %-10s %-12s %-10s %-12s\n', 'p', '|alpha|', 'tau [ns]', 'nu [Hz]', '|nu/numax|');
for p = 1:P
    fprintf('%-4d %-10.4f %-12.1f %-10.2f %-12.3f\n', ...
        p, abs(alphas(p)), taus(p)*1e9, nus(p), nus(p)/nu_max);
end
fprintf('\n');

fprintf('FIM size            : %d x %d\n', size(J_OTFS,1), size(J_OTFS,2));
fprintf('cond(J) OTFS        : %.3e\n', cond(J_OTFS));
fprintf('cond(J) OFDM no-ICI : %.3e\n', cond(J_OFDM));
fprintf('||J_OTFS - J_OFDM||_F / ||J_OFDM||_F = %.3e   (Section 3.2 equivalence)\n', ...
        norm(J_OTFS - J_OFDM, 'fro') / norm(J_OFDM, 'fro'));
fprintf('\n');

fprintf('%-4s %-15s %-15s %-15s %-15s\n', 'p', 'CRB_tau OTFS', 'CRB_tau OFDM', 'CRB_nu OTFS', 'CRB_nu OFDM');
for p = 1:P
    i_tau = 4*(p-1)+3;
    i_nu  = 4*(p-1)+4;
    fprintf('%-4d %-15.3e %-15.3e %-15.3e %-15.3e\n', p, ...
        CRB_OTFS(i_tau,i_tau), CRB_OFDM(i_tau,i_tau), ...
        CRB_OTFS(i_nu, i_nu ), CRB_OFDM(i_nu, i_nu ));
end

% ---------- local helper ----------
function C = invFIM_scaled(J)
    d  = sqrt(diag(J));
    Js = J ./ (d * d.');
    Cs = inv(Js);
    C  = Cs ./ (d * d.');
end
