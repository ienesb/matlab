function taps = dd_response_io(g, cr, use_ici, E, tap_tol, mode)
% DD_RESPONSE_IO  Effective delay-Doppler tap list of the CP-OTFS channel,
% derived from the SAME io-layer model that generates the received signal
% (latex/io/cp_otfs_io.tex eq. OTFS_IO_general / OTFS_IO_ICI), for the MP
% detector's "ideal CSI" arm.
%   taps = dd_response_io(g, cr, use_ici, E, tap_tol, mode)
% A unit DD impulse at (0,0) is passed through io_cp_otfs; the response h[l',k']
% is kept on the path delay rows l_p within +-E Doppler bins of round(k_real_p)
% (E = Inf keeps whole rows) and above tap_tol * max|h|.  Returned as 0-based
% shifts taps.dl, taps.dk and gains taps.gain (S x 1), plus taps.h (N x M).
% mode 'origin' (default) probes one impulse; 'per_row' also probes impulses
% on a few other delay rows and reports the shift-variance residual
% taps.residual (nonzero with ICI, where the DD channel is not a pure 2-D
% convolution).  Differs from Raviteja WCNC'18 eq. (13) / the legacy
% build_effective_dd_eva by (i) the twist phase e^{-j2pi k_real l/(NM)},
% (ii) TS (not T) in the Doppler kernel, and (iii) the SIGN of the fractional-
% Doppler leakage phase: this model gives sum_m e^{j2pi m (k_real - k')/M},
% the printed eq. (13) the conjugate argument (its eqs. 11-12 conjugate its
% own eq. 7 + 9); magnitudes coincide.  Generation and detection share this
% model, so the BER is unaffected (test_rx pins the relation).
    if nargin < 4 || isempty(E), E = 10; end
    if nargin < 5 || isempty(tap_tol), tap_tol = 0.02; end
    if nargin < 6 || isempty(mode), mode = 'origin'; end
    N = g.N;  M = g.M;
    x0 = zeros(N, M);  x0(1, 1) = 1;
    h = io_cp_otfs(x0, g, cr, use_ici);
    mask = false(N, M);
    for p = 1:numel(cr.alpha)
        lp = mod(round(cr.tau(p)/g.Ts), N);
        kc = round(cr.nu(p)*M*g.TS);
        if isinf(E)
            mask(lp+1, :) = true;
        else
            cols = mod(kc + (-E:E), M);
            mask(lp+1, cols+1) = true;
        end
    end
    keep = mask & (abs(h) >= tap_tol*max(abs(h(:))));
    [r, c] = find(keep);
    taps.dl = r - 1;  taps.dk = c - 1;  taps.gain = h(keep);
    taps.h = h;  taps.mask = keep;  taps.S = numel(taps.gain);
    taps.energy_kept = sum(abs(h(keep)).^2) / sum(abs(h(:)).^2);
    if strcmp(mode, 'per_row')
        rows = unique([0, 1, floor(N/3), floor(2*N/3)]);
        resid = 0;
        for l0 = rows(2:end)
            x1 = zeros(N, M);  x1(l0+1, 1) = 1;
            h1 = circshift(io_cp_otfs(x1, g, cr, use_ici), [-l0, 0]);
            resid = max(resid, norm(h1(:) - h(:)) / norm(h(:)));
        end
        taps.residual = resid;
    end
end
