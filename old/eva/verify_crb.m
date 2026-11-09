% verify_crb.m
% CRB / Fisher-information verification on the REDUCED EVA channel.
%
% For each of the four cases {no-ICI, ICI} x {single, multi target}, the FIM is
% computed four ways and checked to agree:
%   (1) OFDM   Monte-Carlo average of the realization FIM  getFIM_OFDM_*
%   (2) OFDM   closed-form expected FIM                    getFIM_OFDM_*_expected
%   (3) OTFS   Monte-Carlo average of the realization FIM  getFIM_CP_OTFS_*
%   (4) OTFS   closed-form expected FIM                    getFIM_CP_OTFS_*_expected
%
% The FIM depends only on second-order data moments, which are (I)SFFT-unitary
% invariant, so (1) = (2) = (3) = (4) up to Monte-Carlo error.  Data is i.i.d.
% unit-energy QPSK: in the time-frequency grid for OFDM, in the delay-Doppler
% grid for OTFS.  The closed forms (2) and (4) are analytically identical (the
% CP-OTFS expected functions delegate to the OFDM ones), so their gap is ~0.
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
rng(1);

%% System parameters
deltaf = 15e3;  T = 1/deltaf;  Ts = T;   % subcarrier spacing / symbol duration [s]
sigma2 = 1.3;                            % noise variance
N = 16;  M = 8;                          % grid size (subcarriers x symbols)
nT = 20000;                              % Monte-Carlo trials

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

tolMC  = 2e-2;                          % MC relative-error tolerance
tolAna = 1e-9;                          % analytic (closed-form) tolerance
relerr = @(A,B) norm(A-B,'fro')/max(norm(B,'fro'),eps);

cases = { 'no-ICI single', false, 1
          'ICI    single', true,  1
          'no-ICI multi ', false, 3
          'ICI    multi ', true,  3 };

fprintf('CRB (FIM) verification on reduced EVA  (N=%d, M=%d, %d MC trials)\n\n', N, M, nT);
fprintf('%-14s | OFDM mc=exp | OTFS mc=exp | OFDM-exp=OTFS-exp | status\n', 'case');
fprintf('%s\n', repmat('-',1,72));

allpass = true;
for ci = 1:size(cases,1)
    name = cases{ci,1};  use_ici = cases{ci,2};  P = cases{ci,3};
    idx = 1:P;
    alpha = alphaR(idx);  tau = tauR(idx);  nu = nuR(idx);

    % ---- closed-form expected FIM (OFDM and OTFS) ----
    Jo_e = fim_ofdm_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P);
    Jt_e = fim_otfs_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P);

    % ---- Monte-Carlo average of the realization FIM ----
    Jo_acc = 0;  Jt_acc = 0;
    for t = 1:nT
        X = qpsk(N,M);    % TF-domain symbols for OFDM
        x = qpsk(N,M);    % DD-domain symbols for OTFS
        Jo_acc = Jo_acc + fim_ofdm_real(X,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P);
        Jt_acc = Jt_acc + fim_otfs_real(x,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P);
    end
    Jo_mc = Jo_acc/nT;  Jt_mc = Jt_acc/nT;

    d_o = relerr(Jo_mc, Jo_e);          % OFDM:  MC average vs expectation
    d_t = relerr(Jt_mc, Jt_e);          % OTFS:  MC average vs expectation
    d_x = relerr(Jo_e, Jt_e);           % cross-modulation expected equality (~0)

    pass = (d_o < tolMC) && (d_t < tolMC) && (d_x < tolAna);
    allpass = allpass && pass;
    fprintf('%-14s | %11.2e | %11.2e | %16.2e | %s\n', ...
            name, d_o, d_t, d_x, tern(pass,'PASS','FAIL'));
end
fprintf('%s\n', repmat('-',1,72));
fprintf('Overall: %s\n', tern(allpass,'PASS','FAIL'));

%% ---- dispatch helpers ----
function J = fim_ofdm_real(X,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P)
    if P == 1
        if use_ici, J = getFIM_OFDM_ICI  (X,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       J = getFIM_OFDM_noICI(X,deltaf,Ts,sigma2,alpha,tau,nu); end
    else
        if use_ici, J = getFIM_OFDM_multi_ICI  (X,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       J = getFIM_OFDM_multi_noICI(X,deltaf,Ts,sigma2,alpha,tau,nu); end
    end
end

function J = fim_otfs_real(x,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P)
    if P == 1
        if use_ici, J = getFIM_CP_OTFS_ICI  (x,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       J = getFIM_CP_OTFS_noICI(x,deltaf,Ts,sigma2,alpha,tau,nu); end
    else
        if use_ici, J = getFIM_CP_OTFS_multi_ICI  (x,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       J = getFIM_CP_OTFS_multi_noICI(x,deltaf,Ts,sigma2,alpha,tau,nu); end
    end
end

function J = fim_ofdm_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P)
    if P == 1
        if use_ici, J = getFIM_OFDM_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,nu);
        else,       J = getFIM_OFDM_noICI_expected(N,M,deltaf,Ts,sigma2,alpha); end
    else
        if use_ici, J = getFIM_OFDM_multi_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       J = getFIM_OFDM_multi_noICI_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu); end
    end
end

function J = fim_otfs_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu,use_ici,P)
    if P == 1
        if use_ici, J = getFIM_CP_OTFS_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       J = getFIM_CP_OTFS_noICI_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu); end
    else
        if use_ici, J = getFIM_CP_OTFS_multi_ICI_expected  (N,M,deltaf,Ts,sigma2,alpha,tau,nu);
        else,       J = getFIM_CP_OTFS_multi_noICI_expected(N,M,deltaf,Ts,sigma2,alpha,tau,nu); end
    end
end

function X = qpsk(N,M)
    b = randi([0 3], N, M);
    X = (1/sqrt(2)) * ((2*bitand(b,1)-1) + 1j*(2*floor(b/2)-1));
end

function s = tern(cond,a,b)
    if cond, s = a; else, s = b; end
end
