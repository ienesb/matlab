% Wrapper: runs plot_mcrb_dd.m, exports figures, prints numerical summary.
close all; clear; clc;

plot_mcrb_dd;

% Save figures
figs = findall(0, 'Type', 'figure');
for i = 1:numel(figs)
    fname = sprintf('plot_mcrb_dd_fig%02d.png', numel(figs) - i + 1);
    exportgraphics(figs(i), fname, 'Resolution', 150);
end
fprintf('Saved %d figures.\n', numel(figs));

% Numerical summary: pull data back from .mat files for diagnostics
r_dd = load('eva_ber_mp_mcrb_dd_resultsv2.mat');
r_tf = load('eva_ber_mp_mcrb_resultsv3.mat');

SNR_dbs   = r_dd.SNR_dbs;
P         = r_dd.P;
SNR_p_dBs = r_dd.SNR_p_dBs;
n_pilot   = r_dd.n_pilot;

ber_csi = mean(r_dd.bers_otfs_csi, 2);
ber_ce  = squeeze(mean(r_dd.bers_otfs_ce, 2));

MCRB_CSI_DD = mean(r_dd.MCRBs_OTFS_CSI, 4);
MCRB_CE_DD  = squeeze(mean(r_dd.MCRBs_OTFS_CE,  4));
MCRB_CSI_TF = mean(r_tf.MCRBs_OTFS,    4);
MCRB_CE_TF  = squeeze(mean(r_tf.MCRBs_OTFS_CE, 4));

IDX_AR  = @(p) 4*p - 3;
IDX_AI  = @(p) 4*p - 2;
IDX_NU  = @(p) 4*p - 1;
IDX_TAU = @(p) 4*p;

fprintf('\n=== Setup ===\n');
fprintf('SNR (dB):   '); fprintf('%6.1f', SNR_dbs); fprintf('\n');
fprintf('SNR_p (dB): '); fprintf('%6.1f', SNR_p_dBs); fprintf('\n');
fprintf('Paths P=%d, n_pilot=%d\n', P, n_pilot);

fprintf('\n=== BER ===\n');
fprintf('SNR     CSI         '); for k=1:n_pilot, fprintf('CE@%2ddB     ', SNR_p_dBs(k)); end; fprintf('\n');
for s = 1:numel(SNR_dbs)
    fprintf('%4.1f  %.3e  ', SNR_dbs(s), ber_csi(s));
    for k = 1:n_pilot, fprintf('%.3e  ', ber_ce(s,k)); end
    fprintf('\n');
end

fprintf('\n=== Delay RMSE [s] (sqrt of MCRB diag) — Path 1 ===\n');
fprintf('SNR  CSI(TF)      CSI(DD)      ');
for k=1:n_pilot, fprintf('CE@%2d(TF)    CE@%2d(DD)    ', SNR_p_dBs(k), SNR_p_dBs(k)); end
fprintf('\n');
for s = 1:numel(SNR_dbs)
    fprintf('%4.1f ', SNR_dbs(s));
    fprintf('%.3e   %.3e   ', sqrt(MCRB_CSI_TF(IDX_TAU(1),IDX_TAU(1),s)), ...
                              sqrt(MCRB_CSI_DD(IDX_TAU(1),IDX_TAU(1),s)));
    for k = 1:n_pilot
        fprintf('%.3e   %.3e   ', sqrt(MCRB_CE_TF(IDX_TAU(1),IDX_TAU(1),s,k)), ...
                                  sqrt(MCRB_CE_DD(IDX_TAU(1),IDX_TAU(1),s,k)));
    end
    fprintf('\n');
end

fprintf('\n=== Doppler RMSE [Hz] — Path 1 (CSI: TF vs DD) ===\n');
for s = 1:numel(SNR_dbs)
    fprintf('SNR=%4.1f  TF=%.3e  DD=%.3e\n', SNR_dbs(s), ...
        sqrt(MCRB_CSI_TF(IDX_NU(1),IDX_NU(1),s)), ...
        sqrt(MCRB_CSI_DD(IDX_NU(1),IDX_NU(1),s)));
end

fprintf('\n=== |alpha| RMSE — Path 1 (CSI: TF vs DD) ===\n');
for s = 1:numel(SNR_dbs)
    var_tf = MCRB_CSI_TF(IDX_AR(1),IDX_AR(1),s) + MCRB_CSI_TF(IDX_AI(1),IDX_AI(1),s);
    var_dd = MCRB_CSI_DD(IDX_AR(1),IDX_AR(1),s) + MCRB_CSI_DD(IDX_AI(1),IDX_AI(1),s);
    fprintf('SNR=%4.1f  TF=%.3e  DD=%.3e\n', SNR_dbs(s), sqrt(var_tf), sqrt(var_dd));
end

fprintf('\n=== Per-path delay RMSE @ SNR=%g dB (CSI, DD) ===\n', SNR_dbs(end));
for p = 1:P
    fprintf('Path %d: tau RMSE=%.3e s, nu RMSE=%.3e Hz, |a| RMSE=%.3e\n', p, ...
        sqrt(MCRB_CSI_DD(IDX_TAU(p),IDX_TAU(p),end)), ...
        sqrt(MCRB_CSI_DD(IDX_NU(p), IDX_NU(p), end)), ...
        sqrt(MCRB_CSI_DD(IDX_AR(p), IDX_AR(p), end) + MCRB_CSI_DD(IDX_AI(p),IDX_AI(p),end)));
end
