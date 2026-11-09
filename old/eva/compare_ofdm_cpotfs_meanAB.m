% compare_ofdm_cpotfs_meanAB.m
% Compares the CLOSED-FORM EXPECTED (mean) MCRB component matrices A and B of
% OFDM vs CP-OTFS, i.i.d. unit-energy QPSK with BSC rate p.  No Monte-Carlo.
%
% Theory (cp_otfs_mcrb.tex, (I)SFFT unitarity):
%   E{A}_cpotfs = E{A}_ofdm                      (exact equality)
%   E{B}_cpotfs = E{B}_ofdm - Cov(b)_ofdm + Cov(b)_cpotfs
%               => E{B}_cpotfs - E{B}_ofdm = Cov(b)_cpotfs - Cov(b)_ofdm
% So A matrices must coincide; B matrices differ only through the bias
% covariance, which is identical on the gain block (no-ICI) but differs on the
% kinematic (tau,nu) blocks and, with ICI, everywhere.
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

deltaf = 15e3;  T = 1/deltaf;  Ts = T;
sigma2 = 1.3;  p = 0.1;
relerr = @(X,Y) norm(X-Y,'fro')/max(norm(Y,'fro'),eps);
fprintf('OFDM vs CP-OTFS expected (mean) A,B   (deltaf=%g, sigma2=%g, p=%.2f)\n\n', deltaf, sigma2, p);

cases = { 'noICI single', @getAB_OFDM_noICI_expected,       @getAB_CP_OTFS_noICI_expected,       1
          'ICI   single', @getAB_OFDM_ICI_expected,         @getAB_CP_OTFS_ICI_expected,         1
          'noICI multi ', @getAB_OFDM_multi_noICI_expected, @getAB_CP_OTFS_multi_noICI_expected, 2
          'ICI   multi ', @getAB_OFDM_multi_ICI_expected,   @getAB_CP_OTFS_multi_ICI_expected,   2 };

N = 16; M = 8;
for ci = 1:size(cases,1)
    name = cases{ci,1};  fO = cases{ci,2};  fC = cases{ci,3};  P = cases{ci,4};
    if P == 1
        alpha = 0.8*exp(1j*0.7);  tau = 370e-9;  nu = 300;
    else
        alpha = [1; 0.7*exp(1j*1.1)];  tau = [0; 150e-9];  nu = [300; -200];
    end

    [Ao, Bo] = fO(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p);
    [Ac, Bc] = fC(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p);

    dA = Ac - Ao;  dB = Bc - Bo;
    fprintf('%s : ||A_cp - A_ofdm|| rel = %.2e   |   ||B_cp - B_ofdm|| rel = %.2e\n', ...
            name, relerr(Ac,Ao), relerr(Bc,Bo));

    % Where does B differ?  Block-resolved (first target only).
    lab = {'aR','aI','tau','nu'};
    fprintf('        B difference (|Bc-Bo|, target-1 4x4 block):\n');
    for r = 1:4
        fprintf('          %-3s ', lab{r});
        for c = 1:4
            fprintf('%9.2e ', abs(dB(r,c)));
        end
        fprintf('\n');
    end
    fprintf('        max|A_cp - A_ofdm| = %.2e (should be ~0)\n\n', max(abs(dA(:))));
end
