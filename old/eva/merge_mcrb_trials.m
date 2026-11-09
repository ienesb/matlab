function merge_mcrb_trials(fileA, fileB, fileout)
% MERGE_MCRB_TRIALS  Combine two MCRB runs over the SAME SNR grid / config into a
% single result with the pooled trial count (e.g. 100 + 100 -> 200 trials).  The
% stored G_bar/F_bar/BB_bar/Gc_bar/ber_data are trial-AVERAGES, so the pooled
% average is the trial-count-weighted mean; kappa_max is the max; and the RMSE
% curves are RECOMPUTED from the pooled matrices (RMSE is nonlinear -> cannot be
% averaged directly).
%
%   merge_mcrb_trials                       % _100a + _100b -> mcrb_ce_ici_results.mat
%   merge_mcrb_trials(A, B, out)            % explicit files
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(fileA),   fileA   = fullfile(thisdir,'results','mcrb_ce_ici_results_100a.mat'); end
    if nargin < 2 || isempty(fileB),   fileB   = fullfile(thisdir,'results','mcrb_ce_ici_results_100b.mat'); end
    if nargin < 3 || isempty(fileout), fileout = fullfile(thisdir,'results','mcrb_ce_ici_results.mat'); end

    A = load(fileA);  B = load(fileB);
    assert(isequal(A.SNR_dbs,B.SNR_dbs) && isequal(A.use_full_eva,B.use_full_eva) && ...
           A.P==B.P && isequal(A.SNRp_dB,B.SNRp_dB) && A.N==B.N && A.M==B.M && ...
           isequal(A.configs,B.configs), ...
           'runs must share SNR grid / channel / config to pool trials.');

    nA = A.nTrials;  nB = B.nTrials;  nT = nA + nB;
    wavg = @(a,b) (nA*a + nB*b) / nT;                 % trial-count-weighted mean

    G_bar  = wavg(A.G_bar,  B.G_bar );
    F_bar  = wavg(A.F_bar,  B.F_bar );
    BB_bar = wavg(A.BB_bar, B.BB_bar);
    Gc_bar = wavg(A.Gc_bar, B.Gc_bar);
    ber_data  = wavg(A.ber_data, B.ber_data);
    kappa_max = max(A.kappa_max, B.kappa_max);

    % ---- recompute RMSE curves from the pooled matrices ----
    SNR_dbs=A.SNR_dbs; configs=A.configs; idx_tau=A.idx_tau; idx_nu=A.idx_nu;
    nSNR=numel(SNR_dbs); nC=size(configs,1);
    rmse_tau_mcrb=zeros(nSNR,nC); rmse_nu_mcrb=zeros(nSNR,nC);
    rmse_tau_crb =zeros(nSNR,nC); rmse_nu_crb =zeros(nSNR,nC);
    for ci=1:nC
        for si=1:nSNR
            N0 = 1/db2pow(SNR_dbs(si));
            LB  = N0*G_bar(:,:,si,ci) + F_bar(:,:,si,ci) + BB_bar(:,:,si,ci);
            CRB = N0*Gc_bar(:,:,si,ci);
            rmse_tau_mcrb(si,ci)=sqrt(sum(real(diag(LB(idx_tau,idx_tau)))));
            rmse_nu_mcrb(si,ci) =sqrt(sum(real(diag(LB(idx_nu, idx_nu )))));
            rmse_tau_crb(si,ci) =sqrt(sum(real(diag(CRB(idx_tau,idx_tau)))));
            rmse_nu_crb(si,ci)  =sqrt(sum(real(diag(CRB(idx_nu, idx_nu )))));
        end
    end

    % ---- carry over shared metadata, set pooled trial count ----
    SNRp_dB=A.SNRp_dB; use_full_eva=A.use_full_eva; evaname=A.evaname;
    N=A.N; M=A.M; P=A.P; bin_tau=A.bin_tau; bin_nu=A.bin_nu;  nTrials=nT;

    save(fileout, ...
         'SNR_dbs','SNRp_dB','configs','use_full_eva','evaname','N','M','P','nTrials', ...
         'G_bar','F_bar','BB_bar','Gc_bar','kappa_max','ber_data','idx_tau','idx_nu', ...
         'rmse_tau_mcrb','rmse_nu_mcrb','rmse_tau_crb','rmse_nu_crb','bin_tau','bin_nu');
    fprintf('Pooled %d + %d = %d trials over SNR=[%s] dB -> %s\n', ...
            nA, nB, nT, num2str(SNR_dbs), 'mcrb_ce_ici_results.mat');

    plot_mcrb_ce_ici(fileout);
end
