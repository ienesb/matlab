function [x, masks] = otfs_pilot_frame(Xd, g, ch, x_p, guard_extra)
% OTFS_PILOT_FRAME  Embedded pilot + full-Doppler delay guard (Raviteja TVT'19
% eq. 8, transposed to delay rows): pilot x_p at (l_p, k_p), zeros on ALL Doppler
% bins of the guard rows 0..G-1, data elsewhere.
%   [x, masks] = otfs_pilot_frame(Xd, g, ch, x_p, guard_extra)
% guard_extra (default ch.l_tau) sets G = l_tau + 1 + guard_extra and l_p =
% guard_extra; the default is the full guard G = 2 l_tau + 1, l_p = l_tau
% (5248 cells = 8.0 % at N = 512, M = 128, l_tau = 20).  Shrinking guard_extra
% is the receiver knob studied in old/eva/mcrb_sweep_grid.m.
% masks: data, known (pilot + guard), known_val, ce_rows (1-based rows
% l_p..l_p+l_tau), l_p, k_p, G, n_nondata.
    [N, M] = size(Xd);
    l_tau = ch.l_tau;
    if nargin < 5 || isempty(guard_extra), guard_extra = l_tau; end
    assert(guard_extra >= 0 && l_tau + 1 + guard_extra <= N, 'otfs_pilot_frame:guard', 'guard does not fit');
    G   = l_tau + 1 + guard_extra;
    l_p = guard_extra;
    k_p = floor(M/2);
    x = Xd;
    x(1:G, :) = 0;
    x(l_p+1, k_p+1) = x_p;
    masks.data = true(N, M);  masks.data(1:G, :) = false;
    masks.known = ~masks.data;
    masks.known_val = zeros(N, M);  masks.known_val(l_p+1, k_p+1) = x_p;
    masks.ce_rows = (l_p+1):(l_p+l_tau+1);
    masks.l_p = l_p;  masks.k_p = k_p;  masks.G = G;
    masks.n_nondata = G*M;
    masks.x_p = x_p;
end
