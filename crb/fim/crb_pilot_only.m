function S = crb_pilot_only(g, cr, use_ici, pilot_frac, sigma2, wf)
% CRB_PILOT_ONLY  Genie CRB from a known comb of pilot rows only (no data),
% the "pilot-only" reference curve of the MCRB figures.  A fraction pilot_frac
% (default 0.20) of the N subcarriers carries unit pilots on every symbol.
% No mismatch is possible, so this is a pure CRB that scales as sigma2.
    if nargin < 4 || isempty(pilot_frac), pilot_frac = 0.20; end
    if nargin < 5 || isempty(sigma2), sigma2 = 1; end
    if nargin < 6 || isempty(wf), wf = 'cp-ofdm'; end
    Npil = round(pilot_frac*g.N);
    pilot_idx = unique(round(linspace(1, g.N, Npil)));
    Xpil = zeros(g.N, g.M);  Xpil(pilot_idx, :) = 1;
    J = fim(Xpil, wf, g, cr, use_ici, sigma2);
    S = crb_from_fim(J, g, cr);
    S.pilot_idx = pilot_idx;
end
