function out = diag_cond_eva(nDraws)
% DIAG_COND_EVA  Conditioning gate for the bound channel (plan Phase 0.8).
% Draws the full 9-path EVA and the 3-tap EVA-reduced channel (Rayleigh and
% constant-magnitude gains), builds the CP-OFDM ICI FIM in bin units, and
% reports the condition number of the kinematic block and the CRB validity.
% Decision rule (README): full EVA is the bound channel only if fewer than 5 %
% of draws exceed kappa = 1e10; otherwise EVA-reduced is the bound workhorse
% and full EVA stays the BER channel.  Saves results/diag_cond_eva.mat and
% figures/fig_cond_eva.png.
    if nargin < 1, nDraws = 500; end
    root = fileparts(fileparts(mfilename('fullpath')));
    N = 512; M = 128; df = 15e3;  Ts = 1/(N*df);
    g = make_grid(N, M, df, 20*Ts, 'cp');
    X = gen_symbols(N, M, '4qam', crn_stream(9001, 0, 1));
    variants = {'full-rayleigh', false, 'rayleigh'; 'full-constmag', false, 'constmag'; ...
                'reduced-rayleigh', true, 'rayleigh'; 'reduced-constmag', true, 'constmag'};
    out = struct();
    fig = new_fig([100 100 1000 700]);
    for v = 1:size(variants, 1)
        ch = chan_eva(g, 'reduced', variants{v,2}, 'gain_model', variants{v,3});
        idx = eta_index(ch.P);
        kappa = zeros(nDraws, 1);  kkin = kappa;  valid = false(nDraws, 1);  trtau = kappa;
        D = eta_scale(g, ch.P);
        parfor t = 1:nDraws
            cr = chan_realize(ch, g, crn_stream(9000, v, t));
            J = fim(X, 'cp-ofdm', g, cr, true, 1);
            S = crb_from_fim(J, g, cr);
            kappa(t) = S.kappa;  valid(t) = S.valid;  trtau(t) = S.tr_tau;
            Jb = D*J*D;  kkin(t) = cond(Jb(idx.kin, idx.kin));
        end
        name = strrep(variants{v,1}, '-', '_');
        out.(name) = struct('kappa', kappa, 'kappa_kin', kkin, 'valid', valid, 'tr_tau', trtau, ...
            'drop_frac_1e10', mean(kappa > 1e10), 'drop_frac_1e8', mean(kappa > 1e8), 'median_kappa', median(kappa));
        fprintf('%-18s  median kappa %.2e   P90 %.2e   max %.2e   >1e8: %.1f%%   >1e10: %.1f%%   invalid: %.1f%%\n', ...
            variants{v,1}, median(kappa), prctile(kappa, 90), max(kappa), 100*mean(kappa > 1e8), 100*mean(kappa > 1e10), 100*mean(~valid));
        subplot(2, 2, v);
        histogram(log10(kappa), 30);  xlabel('log_{10} cond(J_{bin})');  ylabel('draws');
        title(sprintf('%s: >1e10 in %.1f%%', variants{v,1}, 100*mean(kappa > 1e10)), 'Interpreter', 'none');
        grid on;
    end
    fr = out.full_rayleigh.drop_frac_1e10;
    if fr < 0.05
        out.decision = 'full EVA (9 paths) is the bound channel';
    else
        out.decision = 'EVA-reduced (3 taps) is the bound workhorse; full EVA stays the BER channel';
    end
    fprintf('DECISION: %s (full-Rayleigh drop fraction at kappa > 1e10: %.1f%%)\n', out.decision, 100*fr);
    out.nDraws = nDraws;  out.rule = 'full EVA if < 5% of Rayleigh draws exceed kappa 1e10';
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    save(fullfile(root, 'results', 'diag_cond_eva.mat'), '-struct', 'out');
    print(fig, fullfile(root, 'figures', 'fig_cond_eva.png'), '-dpng', '-r120');
    close(fig);
end
