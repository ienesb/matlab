function fig = plot_mcrb_vs_p(res, statistic)
% PLOT_MCRB_VS_P  Delay and Doppler RMSE bounds vs BSC rate p at each SNR of a
% run_mcrb_bsc_eva result: CRB (dotted), MCRB0 (solid), LB = MCRB0 + BB (bold),
% sqrt(BB) (thin dash-dot), MCRB@eta_bar (dashed), mean-field reference
% CRB/(1-2p) (thin grey).  The statistic ('mean' | 'median') is taken over
% valid trials of the PER-TRIAL composite bounds (mcrb_stats).
    if nargin < 2, statistic = res.opts.statistic; end
    root = fileparts(fileparts(mfilename('fullpath')));
    p = res.p_list;  nW = numel(res.waveforms);  nI = numel(res.ici_list);  nS = numel(res.snr_db);
    cols = lines(nW*nI);
    fig = new_fig([60 60 420*nS 720]);
    for is = 1:nS
        N0 = 10^(-res.snr_db(is)/10);
        for blk = 1:2
            if blk == 1, tag = 'tau'; sc = 1e9; lab = 'delay RMSE bound (ns)'; else, tag = 'nu'; sc = 1; lab = 'Doppler RMSE bound (Hz)'; end
            S = mcrb_stats(res.(['G_' tag]), res.(['F_' tag]), res.(['BB_' tag]), res.(['Gc_' tag]), res.valid, N0, statistic);
            T = mcrb_stats(res.(['Gt_' tag]), res.(['Ft_' tag]), 0*res.(['BB_' tag]), res.(['Gc_' tag]), res.valid, N0, statistic);
            subplot(2, nS, (blk-1)*nS + is);  hold on;  k = 0;
            for iw = 1:nW, for ii = 1:nI
                k = k + 1;  c = cols(k, :);
                nm = sprintf('%s %s', res.waveforms{iw}, ternary(res.ici_list(ii), 'ICI', 'no-ICI'));
                semilogy(p, sc*squeeze(S.rmse_crb(1, :, iw, ii)),  ':', 'Color', c, 'LineWidth', 1.0, 'DisplayName', [nm ' CRB']);
                semilogy(p, sc*squeeze(S.rmse_mcrb(1, :, iw, ii)), '-', 'Color', c, 'LineWidth', 1.2, 'DisplayName', [nm ' MCRB0']);
                semilogy(p, sc*squeeze(S.rmse_lb(1, :, iw, ii)),   '-', 'Color', c, 'LineWidth', 2.2, 'DisplayName', [nm ' LB']);
                semilogy(p, sc*squeeze(S.rmse_bb(1, :, iw, ii)),   '-.', 'Color', c, 'LineWidth', 0.8, 'DisplayName', [nm ' sqrt(BB)']);
                semilogy(p, sc*squeeze(T.rmse_mcrb(1, :, iw, ii)), '--', 'Color', c, 'LineWidth', 1.0, 'DisplayName', [nm ' MCRB@eta\_bar']);
                if k == 1
                    semilogy(p, sc*squeeze(S.rmse_crb(1, :, iw, ii))./(1 - 2*p), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'DisplayName', 'mean-field CRB/(1-2p)');
                end
            end, end
            set(gca, 'YScale', 'log');  grid on;  xlabel('BSC rate p');  ylabel(lab);
            title(sprintf('SNR = %g dB', res.snr_db(is)));
            if is == 1 && blk == 1, legend('Location', 'best', 'FontSize', 5, 'NumColumns', 2); end
        end
    end
    sgtitle(sprintf('%s: %d trials, %s over valid trials (drop <= %.1f%%), %s gains, Tcp = %d Ts', res.ch.name, res.nT, statistic, ...
        100*max(res.drop_frac(:)), res.ch.gain_model, res.opts.Tcp_samples), 'Interpreter', 'none');
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    print(fig, fullfile(root, 'figures', sprintf('fig_mcrb_vs_p_%s_%s%s_%s.png', lower(strrep(res.ch.name, '-', '_')), res.opts.tier, res.opts.tag, statistic)), '-dpng', '-r120');
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
