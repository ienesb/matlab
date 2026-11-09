% test_getab_eta.m
% Verifies the optional eta_eval evaluation point added to the getAB wrappers
% (OFDM and CP-OTFS), for the four cases {no-ICI, ICI} x {single, multi target}.
%   (1) eta_eval = ground truth  ==  default (no eta_eval)   [backward compat]
%   (2) eta_eval = pseudo-true   ==  getAB_OFDM_pseudoTrue    [OFDM]
%   (3) getAB_CP_OTFS_*(x,xhat,.,eta_eval) == getAB_OFDM_*(ISFFT x, ISFFT xhat,.,eta_eval)
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
rng(7);
deltaf = 15e3;  Ts = 1/deltaf;  sigma2 = 1.3;  p = 0.1;
fc = 2e9; v = 120/3.6; c = 3e8; numax = v*fc/c;
N = 64;  M = 32;
tauR = [0;370;150]*1e-9;  nuR = numax*cos([0.5;2.0;3.5]);
alR  = 10.^([0;-0.6;-1.4]/20).*exp(1j*[0.7;1.1;-0.5]);
cases = {'no-ICI single',false,1;'ICI    single',true,1;'no-ICI multi ',false,3;'ICI    multi ',true,3};
relerr = @(a,b) norm(a-b,'fro')/max(norm(b,'fro'),eps);
packEta = @(al,ta,nv) reshape([real(al(:)).'; imag(al(:)).'; ta(:).'; nv(:).'], [], 1);

fprintf('%-14s | eta=gt vs default | OFDM eta0 vs pseudoTrue | OTFS vs OFDM(ISFFT)\n','case');
for ci = 1:4
    use_ici = cases{ci,2};  P = cases{ci,3};  idx = 1:P;
    alpha = alR(idx);  tau = tauR(idx);  nu = nuR(idx);
    eta_gt = packEta(alpha, tau, nu);
    X = qpsk(N,M);  Xh = bsc_qpsk(X,p);
    x = qpsk(N,M);  xh = bsc_qpsk(x,p);

    if P == 1
        if use_ici, fO=@getAB_OFDM_ICI;        fC=@getAB_CP_OTFS_ICI;
        else,       fO=@getAB_OFDM_noICI;      fC=@getAB_CP_OTFS_noICI; end
    else
        if use_ici, fO=@getAB_OFDM_multi_ICI;  fC=@getAB_CP_OTFS_multi_ICI;
        else,       fO=@getAB_OFDM_multi_noICI; fC=@getAB_CP_OTFS_multi_noICI; end
    end

    % (1) default vs explicit ground-truth eta_eval
    [A0,B0] = fO(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu);
    [A1,B1] = fO(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,eta_gt);
    d_gt = max(relerr(A0,A1), relerr(B0,B1));

    % (2) eta_eval = pseudo-true  vs  getAB_OFDM_pseudoTrue
    eta0 = pseudoTrue_OFDM(X,Xh,deltaf,Ts,alpha,tau,nu,use_ici);
    [Ap,Bp]   = fO(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,eta0);
    [Apt,Bpt] = getAB_OFDM_pseudoTrue(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici);
    d_pt = max(relerr(Ap,Apt), relerr(Bp,Bpt));

    % (3) CP-OTFS wrapper with eta_eval  vs  OFDM at ISFFT-mapped grids
    eta0o = pseudoTrue_CP_OTFS(x,xh,deltaf,Ts,alpha,tau,nu,use_ici);
    [Ao,Bo]   = fC(x,xh,deltaf,Ts,sigma2,alpha,tau,nu,eta0o);
    [Aof,Bof] = fO(isfft_u(x),isfft_u(xh),deltaf,Ts,sigma2,alpha,tau,nu,eta0o);
    d_otfs = max(relerr(Ao,Aof), relerr(Bo,Bof));

    fprintf('%-14s | %14.1e | %18.1e | %15.1e\n', cases{ci,1}, d_gt, d_pt, d_otfs);
end

function X = qpsk(N,M)
    b = randi([0 3], N, M);  X = (1/sqrt(2))*((2*bitand(b,1)-1)+1j*(2*floor(b/2)-1));
end
function Xh = bsc_qpsk(X,p)
    sR=sign(real(X)); sI=sign(imag(X));
    fR=1-2*(rand(size(X))<p); fI=1-2*(rand(size(X))<p);
    Xh=(sR.*fR+1j*(sI.*fI))/sqrt(2);
end
