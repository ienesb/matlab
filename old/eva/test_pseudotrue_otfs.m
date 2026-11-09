% test_pseudotrue_otfs.m
% Sanity-checks pseudoTrue_CP_OTFS for the four cases {no-ICI, ICI} x {single,
% multi target} on the reduced EVA channel.
%   (1) Matched xhat = x  -> eta0 = eta_bar exactly (resnorm ~ 0).
%   (2) Consistency: pseudoTrue_CP_OTFS(x,xhat) == pseudoTrue_OFDM(ISFFT x, ISFFT xhat)
%       (the unitarity delegation, machine precision).
%   (3) BSC (delay-Doppler) -> tau0 ~ tau_bar, |alpha0| < |alpha_bar|.
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
rng(3);
deltaf = 15e3;  Ts = 1/deltaf;  p = 0.1;
fc = 2e9; v = 120/3.6; c = 3e8; numax = v*fc/c;
N = 64;  M = 32;
tauR = [0;370;150]*1e-9;  nuR = numax*cos([0.5;2.0;3.5]);
alR  = 10.^([0;-0.6;-1.4]/20).*exp(1j*[0.7;1.1;-0.5]);
cases = {'no-ICI single',false,1;'ICI    single',true,1;'no-ICI multi ',false,3;'ICI    multi ',true,3};
relerr = @(a,b) norm(a-b)/max(norm(b),eps);

for ci = 1:4
    use_ici = cases{ci,2};  P = cases{ci,3};  idx = 1:P;
    alpha = alR(idx);  tau = tauR(idx);  nu = nuR(idx);
    x = qpsk(N,M);

    % (1) matched
    [~,~,t0,~,rn] = pseudoTrue_CP_OTFS(x, x, deltaf, Ts, alpha, tau, nu, use_ici);
    matched = relerr([t0], [tau]) + rn;

    % (2) consistency with OFDM at ISFFT-mapped grids
    xh = bsc_qpsk(x, p);
    [e_otfs] = pseudoTrue_CP_OTFS(x, xh, deltaf, Ts, alpha, tau, nu, use_ici);
    [e_ofdm] = pseudoTrue_OFDM(isfft_u(x), isfft_u(xh), deltaf, Ts, alpha, tau, nu, use_ici);
    consistency = relerr(e_otfs, e_ofdm);

    % (3) BSC gain attenuation / delay shift
    [~, a0, t0] = pseudoTrue_CP_OTFS(x, xh, deltaf, Ts, alpha, tau, nu, use_ici);

    fprintf('=== %s ===\n', cases{ci,1});
    fprintf('  matched (resnorm+rel.err): %.2e\n', matched);
    fprintf('  consistency vs OFDM(ISFFT): %.2e\n', consistency);
    fprintf('  BSC: |a0(1)|=%.3f (gt %.3f, (1-2p)gt %.3f)  dtau1=%+.2f ns\n', ...
            abs(a0(1)), abs(alpha(1)), (1-2*p)*abs(alpha(1)), (t0(1)-tau(1))/1e-9);
end

function X = qpsk(N,M)
    b = randi([0 3], N, M);  X = (1/sqrt(2))*((2*bitand(b,1)-1)+1j*(2*floor(b/2)-1));
end
function Xh = bsc_qpsk(X,p)
    sR=sign(real(X)); sI=sign(imag(X));
    fR=1-2*(rand(size(X))<p); fI=1-2*(rand(size(X))<p);
    Xh=(sR.*fR+1j*(sI.*fI))/sqrt(2);
end
