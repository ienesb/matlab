function r = test_mcrb_rx_stats()
% Different solver failures in the two arms must not leak into either
% headline statistic or a paired comparison. Includes legacy MAT metadata.
    r = [];
    o = struct('statistic','mean','kappa_thresh',1e10);
    s = struct('snr_db',0,'waveforms',{{'cp-ofdm'}},'opts',o, ...
        'trn',{{'tr_tau_lb','tr_tau_lb_bsc'}});
    s.trial = struct('tr',[1 4;9 1e8;16 16;25 1e8], ...
        'pcvalid',[1;1;1;1],'exitflag',[1;1;0;1],'kappa',ones(4,1), ...
        'pcvalid_bsc',[1;1;1;0],'exitflag_bsc',[1;0;1;1],'kappa_bsc',ones(4,1), ...
        'nerr',[1;2;3;4],'nbits',100*ones(4,1),'ser',[.01;.02;.03;.04]);
    a = reduce_mcrb_rx_stats(s);
    r = tc(r,'real-arm result excludes only real failures',abs(a.rmse_tau_lb-sqrt(35/3)),1e-14);
    r = tc(r,'control excludes its own failures despite valid real arm',abs(a.rmse_tau_lb_bsc-sqrt(10)),1e-14);
    r = tc(r,'paired comparison uses identical accepted draws',norm([a.stat_mean_paired.rmse_tau_lb a.stat_mean_paired.rmse_tau_lb_bsc]-[1 2]),0);
    r = tc(r,'arm and paired sample counts',norm([a.n_valid a.n_valid_bsc a.n_valid_paired]-[3 2 1]),0);
    r = tc(r,'median applies before square root',abs(a.stat_median.rmse_tau_lb-3),0);
    legacy = s; legacy.trial = rmfield(legacy.trial,'pcvalid_bsc');
    legacy.trial.kappa_bsc(4) = 1e12;
    b = reduce_mcrb_rx_stats(legacy);
    r = tc(r,'legacy control validity reconstructed from saved diagnostics',abs(b.rmse_tau_lb_bsc-sqrt(10)),1e-14);
    legacy.trial.exitflag_bsc(:) = 0; b = reduce_mcrb_rx_stats(legacy);
    r = tc(r,'empty control arm stays missing, never zero',double(~isnan(b.rmse_tau_lb_bsc) || b.n_valid_bsc ~= 0),0);
    s.trial.tr(3,2) = NaN; a = reduce_mcrb_rx_stats(s);
    r = tc(r,'nonfinite control trace excluded from paired population',double(a.n_valid_bsc ~= 1),0);
    s.trial.tr(1,1) = NaN; a = reduce_mcrb_rx_stats(s);
    r = tc(r,'nonfinite real trace invalidates its paired draw',double(a.n_valid ~= 2 || a.n_valid_paired ~= 0),0);
end
