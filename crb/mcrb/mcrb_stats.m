function S = mcrb_stats(G, F, BB, Gc, valid, N0, stat)
% MCRB_STATS  Trial statistics of the PER-TRIAL composite bounds (never of the
% affine pieces separately: the median is not additive and G, F, BB of one trial
% share the same A^-1).  Inputs are per-trial traces (nT x ...) at sigma2 = 1:
%   MCRB0_t = N0 G_t + F_t,   LB_t = MCRB0_t + BB_t,   CRB_t = N0 Gc_t
% valid: logical nT x ... mask (same size);  stat: 'mean' | 'median'.
% Output S (each 1 x ... over the non-trial dimensions): mcrb, lb, bb, crb
% (traces) and rmse_mcrb, rmse_lb, rmse_bb, rmse_crb (sqrt of the clamped
% statistic), plus n_valid.
    sz = size(G);
    mcrb = N0*G + F;  lb = mcrb + BB;  crb = N0*Gc;
    S.mcrb = red(mcrb, valid, stat);  S.lb = red(lb, valid, stat);
    S.bb = red(BB, valid, stat);      S.crb = red(crb, valid, stat);
    S.n_valid = reshape(sum(valid, 1), [1 sz(2:end)]);
    for f = {'mcrb', 'lb', 'bb', 'crb'}
        S.(['rmse_' f{1}]) = sqrt(max(S.(f{1}), 0));
    end
end

function r = red(A, valid, stat)
    sz = size(A);  nT = sz(1);
    A2 = reshape(A, nT, []);  v2 = reshape(valid, nT, []);
    r = nan(1, size(A2, 2));
    for k = 1:size(A2, 2)
        x = A2(v2(:, k), k);  x = x(isfinite(x));
        if isempty(x), continue; end
        if strcmpi(stat, 'median'), r(k) = median(x); else, r(k) = mean(x); end
    end
    r = reshape(r, [1 sz(2:end)]);
end
