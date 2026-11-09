function r = test_regression_old()
% TEST_REGRESSION_OLD  Legacy end-to-end regression (plan 6.3): the old/eva
% mcrb_ce_ici configuration (EVA-reduced, constant-magnitude gains, Tcp = 0,
% dftdd CE, posterior-damped MP 0.3 x 30, 50 trials, trial MEAN) through
% run_mcrb_rx_eva, against old/eva/results/mcrb_ce_ici_results_2000.mat:
% genie delay CRB within 5 %, OFDM BER within 30 %, delay LB within 35 % (the
% trial MEAN of LB is heavy-tailed through BB: the old 50-trial and 2000-trial
% files themselves differ by up to 30 % at 0-12 dB), OTFS BER within 50 % (the
% io-layer channel + io taps give 1.3-1.5x the old DD-circular-model OTFS BER
% at 10-12 dB, a documented model difference) at the shared SNR points.
% Slow (~20 min): part of run_all_tests('slow') only.
    r = [];
    ref.snr = [0 10 12 15 20 25];
    ref.ofdm_lb  = [5.825e-09 4.513e-10 2.8e-10 1.369e-10 6.698e-11 3.634e-11];
    ref.otfs_lb  = [3.519e-09 2.354e-10 1.64e-10 1.145e-10 6.48e-11 3.598e-11];
    ref.ofdm_crb = [6.445e-10 2.043e-10 1.623e-10 1.145e-10 6.48e-11 3.597e-11];
    ref.ofdm_ber = [0.229 0.03465 0.0211 0.01007 0.003308 0.001326];
    ref.otfs_ber = [0.2116 0.003747 0.000282 1.419e-06 2.321e-08 0];
    res = run_mcrb_rx_eva('tier', 50, 'legacy', true, 'snr_db', ref.snr, 'do_bsc', false, 'tag', '_regression');
    io = find(strcmp(res.waveforms, 'cp-ofdm'));  it = find(strcmp(res.waveforms, 'cp-otfs'));
    for k = 1:numel(ref.snr)
        r = tc(r, sprintf('legacy OFDM delay LB at %d dB vs old (35%%)', ref.snr(k)), abs(res.stat_mean.rmse_tau_lb(k, io)/ref.ofdm_lb(k) - 1), 0.35);
        r = tc(r, sprintf('legacy OFDM delay CRB at %d dB vs old (5%%)', ref.snr(k)), abs(res.stat_mean.rmse_tau_crb(k, io)/ref.ofdm_crb(k) - 1), 0.05);
        r = tc(r, sprintf('legacy OTFS delay LB at %d dB vs old (35%%)', ref.snr(k)), abs(res.stat_mean.rmse_tau_lb(k, it)/ref.otfs_lb(k) - 1), 0.35);
        if ref.ofdm_ber(k) > 1e-3, r = tc(r, sprintf('legacy OFDM BER at %d dB vs old (30%%)', ref.snr(k)), abs(res.ber(k, io)/ref.ofdm_ber(k) - 1), 0.30); end
        if ref.otfs_ber(k) > 1e-3, r = tc(r, sprintf('legacy OTFS BER at %d dB vs old (50%%, model difference)', ref.snr(k)), abs(res.ber(k, it)/ref.otfs_ber(k) - 1), 0.50); end
    end
end
