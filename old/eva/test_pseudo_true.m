% test_pseudo_true.m
% Sanity-checks pseudoTrue_OFDM (pseudo_true.tex projection) for the four cases
% {no-ICI, ICI} x {single, multi target} on the reduced EVA channel.
%   (1) Matched data Xhat = X  ->  eta0 = eta_bar exactly (resnorm ~ 0).
%   (2) BSC-corrupted Xhat     ->  tau0 ~ tau_bar, nu0 ~ nu_bar, |alpha0| < |alpha_bar|
%       (gain attenuates toward (1-2p)*alpha_bar; delay/Doppler stay unbiased).
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
rng(3);

deltaf = 15e3;  Ts = 1/deltaf;
fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;
N = 64;  M = 32;  p = 0.1;

% reduced EVA (three strongest taps)
tauR = [0; 370; 150]*1e-9;
nuR  = numax*cos([0.5; 2.0; 3.5]);
alR  = 10.^([0; -0.6; -1.4]/20) .* exp(1j*[0.7; 1.1; -0.5]);

cases = { 'no-ICI single', false, 1
          'ICI    single', true,  1
          'no-ICI multi ', false, 3
          'ICI    multi ', true,  3 };

for ci = 1:size(cases,1)
    name = cases{ci,1};  use_ici = cases{ci,2};  P = cases{ci,3};
    idx = 1:P;
    alpha = alR(idx);  tau = tauR(idx);  nu = nuR(idx);
    X = qpsk(N,M);

    % (1) matched -> exact recovery
    [~, a0, t0, n0, rn] = pseudoTrue_OFDM(X, X, deltaf, Ts, alpha, tau, nu, use_ici);
    err_match = norm([a0;t0;n0] - [alpha;tau;nu]) / norm([alpha;tau;nu]);

    % (2) BSC-corrupted detections
    Xh = bsc_qpsk(X, p);
    [~, a0, t0, n0, rn2] = pseudoTrue_OFDM(X, Xh, deltaf, Ts, alpha, tau, nu, use_ici);

    fprintf('=== %s ===\n', name);
    fprintf('  matched : rel.err(eta0,eta_bar) = %.2e , resnorm = %.2e\n', err_match, rn);
    fprintf('  BSC p=%.2f : resnorm = %.2e\n', p, rn2);
    for q = 1:P
        fprintf(['    tgt %d: |a0|=%.4f (|a_bar|=%.4f, (1-2p)|a_bar|=%.4f) | ' ...
                 'dtau=%+.2f ns  dnu=%+.2f Hz\n'], q, abs(a0(q)), abs(alpha(q)), ...
                 (1-2*p)*abs(alpha(q)), (t0(q)-tau(q))/1e-9, n0(q)-nu(q));
    end
end

function X = qpsk(N,M)
    b = randi([0 3], N, M);
    X = (1/sqrt(2)) * ((2*bitand(b,1)-1) + 1j*(2*floor(b/2)-1));
end
function Xh = bsc_qpsk(X, p)
    sR = sign(real(X)); sI = sign(imag(X));
    fR = 1 - 2*(rand(size(X)) < p);
    fI = 1 - 2*(rand(size(X)) < p);
    Xh = (sR.*fR + 1j*(sI.*fI))/sqrt(2);
end
