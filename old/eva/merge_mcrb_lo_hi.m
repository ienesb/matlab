function merge_mcrb_lo_hi(fileLo, fileHi, fileout)
% MERGE_MCRB_LO_HI  Concatenate two MCRB sweeps over DISJOINT SNR grids (same
% channel/config and trial count) into one full-range result, then (re)draw the
% figure.  Both runs must share use_full_eva, P, nTrials, SNR_p, N, M, configs;
% only SNR_dbs differs.  (For pooling more TRIALS at the same SNR, use
% merge_mcrb_trials instead.)
%
%   merge_mcrb_lo_hi                          % _lo + _hi -> mcrb_ce_ici_results.mat
%   merge_mcrb_lo_hi(A, B, out)              % explicit files
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(fileLo),  fileLo  = fullfile(thisdir,'results','mcrb_ce_ici_results_lo.mat'); end
    if nargin < 2 || isempty(fileHi),  fileHi  = fullfile(thisdir,'results','mcrb_ce_ici_results_hi.mat'); end
    if nargin < 3 || isempty(fileout), fileout = fullfile(thisdir,'results','mcrb_ce_ici_results.mat'); end
    lo = load(fileLo);
    hi = load(fileHi);

    assert(isequal(lo.use_full_eva,hi.use_full_eva) && lo.P==hi.P && ...
           lo.nTrials==hi.nTrials && isequal(lo.SNRp_dB,hi.SNRp_dB) && ...
           lo.N==hi.N && lo.M==hi.M && isequal(lo.configs,hi.configs), ...
           'lo/hi runs have mismatched config -- cannot merge.');

    % ---- merge along the SNR axis, sorted ascending ----
    SNR_dbs = [lo.SNR_dbs(:).' , hi.SNR_dbs(:).'];
    [SNR_dbs, ord] = sort(SNR_dbs);

    G_bar  = cat(3, lo.G_bar,  hi.G_bar );  G_bar  = G_bar(:,:,ord,:);
    F_bar  = cat(3, lo.F_bar,  hi.F_bar );  F_bar  = F_bar(:,:,ord,:);
    BB_bar = cat(3, lo.BB_bar, hi.BB_bar);  BB_bar = BB_bar(:,:,ord,:);
    Gc_bar = cat(3, lo.Gc_bar, hi.Gc_bar);  Gc_bar = Gc_bar(:,:,ord,:);

    catrow = @(a,b) [a; b];
    kappa_max     = catrow(lo.kappa_max,    hi.kappa_max);     kappa_max     = kappa_max(ord,:);
    ber_data      = catrow(lo.ber_data,     hi.ber_data);      ber_data      = ber_data(ord,:);
    rmse_tau_mcrb = catrow(lo.rmse_tau_mcrb,hi.rmse_tau_mcrb);  rmse_tau_mcrb = rmse_tau_mcrb(ord,:);
    rmse_nu_mcrb  = catrow(lo.rmse_nu_mcrb, hi.rmse_nu_mcrb);   rmse_nu_mcrb  = rmse_nu_mcrb(ord,:);
    rmse_tau_crb  = catrow(lo.rmse_tau_crb, hi.rmse_tau_crb);   rmse_tau_crb  = rmse_tau_crb(ord,:);
    rmse_nu_crb   = catrow(lo.rmse_nu_crb,  hi.rmse_nu_crb);    rmse_nu_crb   = rmse_nu_crb(ord,:);

    % ---- carry over the shared metadata ----
    SNRp_dB=hi.SNRp_dB; configs=hi.configs; use_full_eva=hi.use_full_eva;
    evaname=hi.evaname; N=hi.N; M=hi.M; P=hi.P; nTrials=hi.nTrials;
    idx_tau=hi.idx_tau; idx_nu=hi.idx_nu; bin_tau=hi.bin_tau; bin_nu=hi.bin_nu;

    save(fileout, ...
         'SNR_dbs','SNRp_dB','configs','use_full_eva','evaname','N','M','P','nTrials', ...
         'G_bar','F_bar','BB_bar','Gc_bar','kappa_max','ber_data','idx_tau','idx_nu', ...
         'rmse_tau_mcrb','rmse_nu_mcrb','rmse_tau_crb','rmse_nu_crb','bin_tau','bin_nu');
    fprintf('Merged SNR = [%s] dB (%d trials) -> %s\n', num2str(SNR_dbs), nTrials, fileout);

    plot_mcrb_ce_ici(fileout);
end
