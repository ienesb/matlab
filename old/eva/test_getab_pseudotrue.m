% test_getab_pseudotrue.m
% Checks getAB_OFDM_pseudoTrue (A,B at the pseudo-true parameter) and that the
% generalized getAB_OFDM_core is backward compatible (eta_bar = eta by default).
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
rng(5);
deltaf = 15e3;  Ts = 1/deltaf;  sigma2 = 1.3;  p = 0.1;
fc = 2e9; v = 120/3.6; c = 3e8; numax = v*fc/c;
N = 64;  M = 32;
tauR = [0;370;150]*1e-9;  nuR = numax*cos([0.5;2;3.5]);
alR  = 10.^([0;-0.6;-1.4]/20).*exp(1j*[0.7;1.1;-0.5]);
cases = {'no-ICI single',false,1;'ICI    single',true,1;'no-ICI multi ',false,3;'ICI    multi ',true,3};
relerr = @(A,B) norm(A-B,'fro')/max(norm(B,'fro'),eps);

for ci = 1:4
    use_ici = cases{ci,2};  P = cases{ci,3};  idx = 1:P;
    alpha = alR(idx);  tau = tauR(idx);  nu = nuR(idx);
    X = qpsk(N,M);

    % (a) backward compat: core with explicit eta_bar=eta == core default (== old behaviour)
    Xh = bsc_qpsk(X,p);
    [Ad,Bd] = getAB_OFDM_core(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici);
    [Ae,Be] = getAB_OFDM_core(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,alpha,tau,nu);
    compat = max(relerr(Ad,Ae), relerr(Bd,Be));

    % (b) matched data -> pseudo-true == truth -> A=-B=-J
    [A0,B0] = getAB_OFDM_pseudoTrue(X,X,deltaf,Ts,sigma2,alpha,tau,nu,use_ici);
    matched = relerr(-A0,B0);

    % (c) BSC: A,B at pseudo-true vs at ground truth
    [Apt,Bpt,~,a0,t0,n0] = getAB_OFDM_pseudoTrue(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici);

    fprintf('=== %s ===\n', cases{ci,1});
    fprintf('  backward-compat (explicit vs default eta_bar): %.2e\n', compat);
    fprintf('  matched ||A+B||/||B|| (=CRB check): %.2e\n', matched);
    fprintf('  BSC: A(pt vs gt)=%.2e  B(pt vs gt)=%.2e | |a0(1)|=%.3f(gt %.3f) dtau1=%+.2fns\n', ...
            relerr(Apt,Ad), relerr(Bpt,Bd), abs(a0(1)), abs(alpha(1)), (t0(1)-tau(1))/1e-9);
end

function X = qpsk(N,M)
    b = randi([0 3], N, M);  X = (1/sqrt(2))*((2*bitand(b,1)-1)+1j*(2*floor(b/2)-1));
end
function Xh = bsc_qpsk(X,p)
    sR=sign(real(X)); sI=sign(imag(X));
    fR=1-2*(rand(size(X))<p); fI=1-2*(rand(size(X))<p);
    Xh=(sR.*fR+1j*(sI.*fI))/sqrt(2);
end
