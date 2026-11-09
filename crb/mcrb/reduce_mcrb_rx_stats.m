function res = reduce_mcrb_rx_stats(res, o)
% REDUCE_MCRB_RX_STATS  Validate each detector/control arm independently.
% Real-only summaries retain the real-arm population. Randomized-control
% summaries use its own validity. stat_{mean,median}_paired uses the same
% intersection for every trace, for real/control comparisons and uncertainty.
% Historical MAT files without pcvalid_bsc use the saved control exit flag,
% the same intrinsic 1e12 condition gate, and finite control traces; original
% files need not be rewritten to recompute their summaries.
    if nargin < 2, o = res.opts; end
    nS = numel(res.snr_db); nW = numel(res.waveforms); trn = res.trn; nQ = numel(trn);
    T = res.trial;
    iscontrol = endsWith(trn, '_bsc');
    finite_real = reshape(all(isfinite(T.tr(:, ~iscontrol & ~endsWith(trn,'_pil'), :, :)), 2), size(T.kappa));
    valid = T.pcvalid > 0 & T.exitflag > 0 & isfinite(T.kappa) & T.kappa <= o.kappa_thresh & finite_real;
    finite_control = reshape(all(isfinite(T.tr(:, iscontrol, :, :)), 2), size(T.kappa));
    if isfield(T, 'pcvalid_bsc')
        control_pcvalid = T.pcvalid_bsc > 0;
        res.control_validity_source = 'saved pcvalid_bsc';
    else
        control_pcvalid = T.kappa_bsc < 1e12 & finite_control;
        res.control_validity_source = 'legacy reconstruction: finite control traces and kappa_bsc < 1e12';
    end
    valid_bsc = control_pcvalid & T.exitflag_bsc > 0 & isfinite(T.kappa_bsc) & T.kappa_bsc <= o.kappa_thresh & finite_control;
    valid_pair = valid & valid_bsc;
    res.trial.valid = valid; res.trial.valid_bsc = valid_bsc; res.trial.valid_paired = valid_pair;
    for f = {'drop_frac','drop_frac_bsc','drop_frac_paired','kappa_max','ber','ser','n_valid','n_valid_bsc','n_valid_paired'}
        res.(f{1}) = nan(nS, nW);
    end
    for st = {'mean','median'}
        for q = 1:nQ
            fn = strrep(trn{q}, 'tr_', 'rmse_');
            res.(['stat_' st{1}]).(fn) = nan(nS,nW);
            res.(['stat_' st{1} '_paired']).(fn) = nan(nS,nW);
        end
    end
    for iw = 1:nW, for si = 1:nS
        if ~any(isfinite(T.tr(:,1,si,iw))), continue; end
        vr = valid(:,si,iw); vb = valid_bsc(:,si,iw); vp = valid_pair(:,si,iw);
        res.drop_frac(si,iw) = 1-mean(vr); res.drop_frac_bsc(si,iw) = 1-mean(vb); res.drop_frac_paired(si,iw) = 1-mean(vp);
        res.n_valid(si,iw) = nnz(vr); res.n_valid_bsc(si,iw) = nnz(vb); res.n_valid_paired(si,iw) = nnz(vp);
        res.kappa_max(si,iw) = max(T.kappa(:,si,iw));
        res.ber(si,iw) = sum(T.nerr(:,si,iw))/max(sum(T.nbits(:,si,iw)),1);
        res.ser(si,iw) = mean(T.ser(:,si,iw));
        for q = 1:nQ
            v = vr; if iscontrol(q), v = vb; end
            fn = strrep(trn{q}, 'tr_', 'rmse_');
            for st = {'mean','median'}
                res.(['stat_' st{1}]).(fn)(si,iw) = root_stat(T.tr(v,q,si,iw), st{1});
                res.(['stat_' st{1} '_paired']).(fn)(si,iw) = root_stat(T.tr(vp,q,si,iw), st{1});
            end
        end
    end, end
    hs = res.(['stat_' o.statistic]); fn = fieldnames(hs);
    for k = 1:numel(fn), res.(fn{k}) = hs.(fn{k}); end
    res.statistic = sprintf('%s of per-trial bound traces over each arm''s valid trials; paired summaries use their intersection; RMSE = sqrt(statistic of trace), kappa <= %g', o.statistic, o.kappa_thresh);
end

function y = root_stat(x, st)
    x = x(isfinite(x));
    if isempty(x), y = NaN; return; end
    if strcmp(st,'median'), v = median(x); else, v = mean(x); end
    y = sqrt(max(v,0));
end
