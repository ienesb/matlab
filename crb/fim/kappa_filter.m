function out = kappa_filter(stats, mode, thresh, p)
% KAPPA_FILTER  Omega_N-type trial filters; always reports the drop fraction.
%   out = kappa_filter(stats, mode, thresh, p)
% mode 'cond'  : stats.kappa (nT x 1) condition numbers of the bin-unit A (or J);
%                valid = kappa <= thresh   (thresh: 1e4 MT40, Inf EVA)
% mode 'liu'   : stats.J (nT cell of matrices), stats.Jbar; Liu Prop. 2
%                rho_i = ||Jbar^{-1} (J_i - Jbar)||_2;  valid = rho_i < thresh (1)
% mode 'kv'    : kappa_validity.tex Sec. "Measuring kappa directly":
%                stats.Akin (nT cell of 2x2 kinematic blocks of A),
%                stats.Jkin (2x2 kinematic block of the genie Jbar);
%                Abar = mean, kappa_i = ||Jkin^{-1}(A_i - Abar)||_2,
%                rho_i = kappa_i/(1-2p); valid = rho_i < thresh (1);
%                also returns kappa_rms, kappa_p90 and the wall p* = (1-kappa_rms)/2.
% out: valid (logical nT x 1), drop_frac, rho (when defined), kappa_rms, kappa_p90, pstar.
    if nargin < 3 || isempty(thresh), thresh = Inf; end
    switch lower(mode)
        case 'cond'
            k = stats.kappa(:);
            out.valid = isfinite(k) & (k <= thresh);
            out.rho = k;
        case 'liu'
            nT = numel(stats.J);  rho = zeros(nT, 1);
            for i = 1:nT, rho(i) = norm(stats.Jbar \ (stats.J{i} - stats.Jbar)); end
            out.valid = rho < thresh;  out.rho = rho;
        case 'kv'
            nT = numel(stats.Akin);
            Abar = zeros(2);
            for i = 1:nT, Abar = Abar + stats.Akin{i}/nT; end
            kap = zeros(nT, 1);
            for i = 1:nT, kap(i) = norm(stats.Jkin \ (stats.Akin{i} - Abar)); end
            rho = kap / (1 - 2*p);
            out.valid = rho < thresh;  out.rho = rho;
            out.kappa_rms = sqrt(mean(kap.^2));
            out.kappa_p90 = prctile_local(kap, 90);
            out.pstar = 0.5*(1 - out.kappa_rms);
        otherwise
            error('kappa_filter:mode', 'unknown mode %s', mode);
    end
    out.drop_frac = 1 - mean(out.valid);
end

function q = prctile_local(x, pct)
    x = sort(x(:));  n = numel(x);
    r = pct/100*(n-1) + 1;  lo = floor(r);  hi = ceil(r);
    q = x(lo) + (r - lo)*(x(hi) - x(lo));
end
