function taps = otfs_threshold_ce(Y, masks, x_p, N0, thr)
% OTFS_THRESHOLD_CE  Embedded-pilot threshold channel estimation, Raviteja
% TVT'19 Sec. III-B1 / eq. (9) (fractional Doppler, full guard): on the pilot
% delay band rows l_p..l_p+l_tau (all Doppler bins), keep |y| >= thr*sqrt(N0)
% and estimate the effective tap y/x_p; map to 0-based shifts relative to the
% pilot cell.  thr = 3 (paper's 3 sigma, equivalently 3 sigma_p on y/x_p).
    if nargin < 5 || isempty(thr), thr = 3; end
    [N, M] = size(Y);
    h = zeros(N, M);
    rows = masks.ce_rows;
    h(rows, :) = Y(rows, :) / x_p;
    h(abs(h) < thr*sqrt(N0)/x_p) = 0;
    [r, c, gc] = find(h);
    taps.dl = mod(r - 1 - masks.l_p, N);
    taps.dk = mod(c - 1 - masks.k_p, M);
    taps.gain = gc;
    taps.S = numel(gc);
    taps.h = h;
end
