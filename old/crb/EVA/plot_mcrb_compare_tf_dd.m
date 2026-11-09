% plot_mcrb_compare_tf_dd.m
% Plot delay-MCRB comparison between OTFS (getMCRB_DD) and OFDM (getMCRB_TF)
% from eva_mcrb_compare_tf_dd.m output.  One figure per path.
clear;

r1 = load('eva_mcrb_fixed_ber_resultsv5.mat');
r2 = load('eva_mcrb_fixed_ber_resultsv5.mat');
r3 = load('eva_mcrb_fixed_ber_resultsv5.mat');

SNR_dbs     = r1.SNR_dbs;
P           = r1.P;
target_bers = r1.target_bers;
n_ber       = r1.n_ber;

shape = size(r1.MCRBs_OFDM);
s1 = size(r1.MCRBs_OFDM, 4);
s2 = size(r2.MCRBs_OFDM, 4);
s3 = size(r3.MCRBs_OFDM, 4);
shape(4) = s1 + s2 + s3;

% Reduce across MC dimension (median is robust to heavy tails)
MCRB_OTFS = zeros(shape);
MCRB_OTFS(:, :, :, 1:s1, :) = r1.MCRBs_OTFS;
MCRB_OTFS(:, :, :, (s1+1):(s1+s2), :) = r2.MCRBs_OTFS;
MCRB_OTFS(:, :, :, (s1+s2+1):end, :) = r3.MCRBs_OTFS;
MCRB_OTFS = (mean(MCRB_OTFS, 4));   % [4P x 4P x nSNR x 1 x n_ber]

MCRB_OFDM = zeros(shape);
MCRB_OFDM(:, :, :, 1:s1, :) = r1.MCRBs_OFDM;
MCRB_OFDM(:, :, :, (s1+1):(s1+s2), :) = r2.MCRBs_OFDM;
MCRB_OFDM(:, :, :, (s1+s2+1):end, :) = r3.MCRBs_OFDM;
MCRB_OFDM = (mean(MCRB_OFDM, 4));   % [4P x 4P x nSNR x 1 x n_ber]

colors  = lines(n_ber);
ber_lab = arrayfun(@(b) sprintf('BER=%.3g', b), target_bers, 'UniformOutput', false);

% Delay parameter index for path p:  eta = [a_R, a_I, nu, tau] -> tau = 4p
IDX_TAU = @(p) 4*p;

for p = 1:P
    figure;
    rmse_otfs = sqrt(squeeze(MCRB_OTFS(IDX_TAU(p), IDX_TAU(p), 1, 1, :)));
    rmse_ofdm = sqrt(squeeze(MCRB_OFDM(IDX_TAU(p), IDX_TAU(p), 1, 1, :)));
    semilogy(target_bers, rmse_otfs,  ...
        'Color', colors(1, :), 'LineWidth', 2); hold on;
    semilogy(target_bers, rmse_ofdm, ...
        'Color', colors(2, :), 'LineWidth', 2);
    grid on;
    xlabel('BER');
    ylabel('Delay RMSE [s]');
    title(sprintf('Path %d — Delay MCRB: OTFS (solid) vs OFDM (dashed)', p));
    legend("OTFS", "OFDM");
end
