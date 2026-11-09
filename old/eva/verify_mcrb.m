% verify_mcrb.m
% MCRB A/B-matrix verification on the REDUCED EVA channel.
%
% For each of the four cases {no-ICI, ICI} x {single, multi target}, the MCRB
% component matrices A and B are computed four ways:
%   (1) OFDM   Monte-Carlo average of the realization  getAB_OFDM_*
%   (2) OFDM   closed-form expectation                 getAB_OFDM_*_expected
%   (3) OTFS   Monte-Carlo average of the realization  getAB_CP_OTFS_*
%   (4) OTFS   closed-form expectation                 getAB_CP_OTFS_*_expected
% where the expected B includes the bias covariance Cov(b):
%   E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b).
%
% Data model: x = i.i.d. unit-energy QPSK, xhat = symbol decisions from a BSC
% with crossover p, applied in the TF grid for OFDM and the DD grid for OTFS.
%
% Equalities checked:
%   A : all four agree -- (1)=(2)=(3)=(4).  E{A} is (I)SFFT-unitary invariant,
%       so the OTFS expected A equals the OFDM expected A exactly.
%   B : MC average = expectation WITHIN each modulation -- (1)=(2) and (3)=(4).
%       Across modulations E{B}_OFDM != E{B}_OTFS in general, because Cov(b)
%       differs by design (cp_otfs_mcrb.tex Sec. 4): the (I)SFFT kernels are not
%       diagonal, so the OTFS Cov(b) does not collapse to a multiple of J.  The
%       cross-modulation B gap is therefore reported for information only.
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
rng(1);

%% System parameters
deltaf = 15e3;  T = 1/deltaf;  Ts = T;   % subcarrier spacing / symbol duration [s]
sigma2 = 1.3;                            % noise variance
p = 0.1;                                 % BSC crossover (symbol-decision error)
N = 16;  M = 8;                          % grid size (subcarriers x symbols)
nT = 40000;                              % Monte-Carlo trials

%% Reduced EVA channel (eva_channel.tex, Sec. "Reduced EVA Channel")
% Three strongest taps of the EVA profile (original indices {0,4,2}).
fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;   % max Doppler ~222.2 Hz
tauR  = [0;   370;  150]*1e-9;          % excess delays [s]
PdB   = [0;  -0.6; -1.4];               % relative powers [dB]
gain  = 10.^(PdB/20);                   % linear tap amplitudes
theta = [0.5; 2.0;  3.5];               % fixed Jakes AoAs -> deterministic Doppler
nuR   = numax*cos(theta);               % per-tap Doppler [Hz]
phiR  = [0.7; 1.1; -0.5];               % fixed gain phases
alphaR = gain .* exp(1j*phiR);          % complex tap gains

tolMC  = 3e-2;                          % MC relative-error tolerance
tolAna = 1e-9;                          % analytic (closed-form) tolerance
relerr = @(A,B) norm(A-B,'fro')/max(norm(B,'fro'),eps);

cases = { 'no-ICI single', false, 1
          'ICI    single', true,  1
          'no-ICI multi ', false, 3
          'ICI    multi ', true,  3 };

fprintf('MCRB (A,B) verification on reduced EVA  (N=%d, M=%d, p=%.2f, %d MC trials)\n\n', N, M, p, nT);
fprintf('%-14s |  A: OFDM mc=exp  OTFS mc=exp  OFDM=OTFS exp |  B: OFDM mc=exp  OTFS mc=exp | (B OFDM-vs-OTFS exp) | status\n', 'case');
fprintf('%s\n', repmat('-',1,118));

allpass = true;
for ci = 1:size(cases,1)
    name = cases{ci,1};  use_ici = cases{ci,2};  P = cases{ci,3};
    idx = 1:P;
    alpha = alphaR(idx);  tau = tauR(idx);  nu = nuR(idx);

    % ---- closed-form expectation (A and B, B includes Cov(b)) ----
    [Ao_e, Bo_e] = ab_ofdm_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p,use_ici,P);
    [At_e, Bt_e] = ab_otfs_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p,use_ici,P);

    % ---- Monte-Carlo average of the realization A and B ----
    Ao_acc=0; Bo_acc=0; At_acc=0; Bt_acc=0;
    for t = 1:nT
        X = qpsk(N,M);   Xh = bsc_qpsk(X,p);    % TF-domain data for OFDM
        x = qpsk(N,M);   xh = bsc_qpsk(x,p);    % DD-domain data for OTFS
        [Ar,Br] = ab_ofdm_real(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P);
        Ao_acc=Ao_acc+Ar; Bo_acc=Bo_acc+Br;
        [Ar,Br] = ab_otfs_real(x,xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P);
        At_acc=At_acc+Ar; Bt_acc=Bt_acc+Br;
    end
    Ao_mc=Ao_acc/nT; Bo_mc=Bo_acc/nT; At_mc=At_acc/nT; Bt_mc=Bt_acc/nT;

    % A: full four-way equality
    dA_o = relerr(Ao_mc, Ao_e);          % OFDM mc vs exp
    dA_t = relerr(At_mc, At_e);          % OTFS mc vs exp
    dA_x = relerr(Ao_e , At_e);          % OFDM-exp vs OTFS-exp (~0)
    % B: mc vs expectation within each modulation
    dB_o = relerr(Bo_mc, Bo_e);          % OFDM mc vs exp
    dB_t = relerr(Bt_mc, Bt_e);          % OTFS mc vs exp
    dB_x = relerr(Bo_e , Bt_e);          % cross-modulation (informational: Cov(b) differs)

    pass = (dA_o<tolMC) && (dA_t<tolMC) && (dA_x<tolAna) && (dB_o<tolMC) && (dB_t<tolMC);
    allpass = allpass && pass;
    fprintf('%-14s | %12.2e %12.2e %12.2e | %12.2e %12.2e | %18.2e | %s\n', ...
            name, dA_o, dA_t, dA_x, dB_o, dB_t, dB_x, tern(pass,'PASS','FAIL'));
end
fprintf('%s\n', repmat('-',1,118));
fprintf('Overall: %s\n', tern(allpass,'PASS','FAIL'));
fprintf(['\nNote: the last column (B OFDM-exp vs OTFS-exp) is expected to be nonzero -- ' ...
         'Cov(b) differs\nbetween OFDM and CP-OTFS (cp_otfs_mcrb.tex Sec. 4); it is not part of the pass criterion.\n']);

%% ---- dispatch helpers ----
function [A,B] = ab_ofdm_real(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P)
    if P == 1
        if use_ici, [A,B] = getAB_OFDM_ICI  (X,Xh,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       [A,B] = getAB_OFDM_noICI(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu); end
    else
        if use_ici, [A,B] = getAB_OFDM_multi_ICI  (X,Xh,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       [A,B] = getAB_OFDM_multi_noICI(X,Xh,deltaf,Ts,sigma2,alpha,tau,nu); end
    end
end

function [A,B] = ab_otfs_real(x,xh,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P)
    if P == 1
        if use_ici, [A,B] = getAB_CP_OTFS_ICI  (x,xh,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       [A,B] = getAB_CP_OTFS_noICI(x,xh,deltaf,Ts,sigma2,alpha,tau,nu); end
    else
        if use_ici, [A,B] = getAB_CP_OTFS_multi_ICI  (x,xh,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       [A,B] = getAB_CP_OTFS_multi_noICI(x,xh,deltaf,Ts,sigma2,alpha,tau,nu); end
    end
end

function [A,B] = ab_ofdm_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p,use_ici,P)
    if P == 1
        if use_ici, [A,B] = getAB_OFDM_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,tau,nu,p);
        else,       [A,B] = getAB_OFDM_noICI_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p); end
    else
        if use_ici, [A,B] = getAB_OFDM_multi_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,tau,nu,p);
        else,       [A,B] = getAB_OFDM_multi_noICI_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p); end
    end
end

function [A,B] = ab_otfs_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p,use_ici,P)
    if P == 1
        if use_ici, [A,B] = getAB_CP_OTFS_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,tau,nu,p);
        else,       [A,B] = getAB_CP_OTFS_noICI_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p); end
    else
        if use_ici, [A,B] = getAB_CP_OTFS_multi_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,tau,nu,p);
        else,       [A,B] = getAB_CP_OTFS_multi_noICI_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,p); end
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

function s = tern(cond,a,b)
    if cond, s = a; else, s = b; end
end
