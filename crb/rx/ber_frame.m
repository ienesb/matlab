function out = ber_frame(wf, g, ch, N0, rx, rs)
% BER_FRAME  One transmit-receive frame: channel draw, symbols, pilots, io-layer
% channel, noise, receiver, hard decisions.  Returns the native-domain grids
% (X, Xhat, with pilots kept) that Phase 3 feeds to the MCRB.
%
%   out = ber_frame(wf, g, ch, N0, rx, rs)
%   wf   'cp-ofdm' | 'cp-otfs'
%   rx   options struct (defaults below); rs a crn_stream (channel, noise and
%        data are drawn in that order so waveforms called with the same stream
%        see the same channel and noise).
% rx fields:
%   modtype 'gray-4qam' etc. ('4qam'), data_amp 1 (symbol amplitude; Raviteja
%   convention: sqrt(SNR_d) with N0 = 1), use_ici true (channel model),
%   pilot 'matched' | 'explicit' | 'none'  (matched: data-referenced budget
%   E_pil = n_nondata * 10^(pilot_boost_dB/10) shared by both waveforms;
%   explicit: x_p given; none: full data frame, ideal CSI only),
%   pilot_boost_dB 0, x_p [], guard_extra [] (= l_tau),
%   ofdm_csi 'ce' | 'ideal1tap' | 'idealbanded', ce_mode 'dftdd', n_dd 2, W_ici 8,
%   otfs_csi 'ce' | 'ideal', thr 3, E 10, tap_tol 0.02, mp_iter 20, mp_damp 0.7,
%   mp_mode 'edge' | 'posterior', mp_tol 1e-4.
% out: bits, bits_hat (data cells), nbits, nerr, ber, ser, X, Xhat, masks, cr,
%      taps (OTFS), H_hat (OFDM), mp_info.
    d = struct('modtype', '4qam', 'data_amp', 1, 'use_ici', true, 'pilot', 'matched', ...
        'pilot_boost_dB', 0, 'x_p', [], 'guard_extra', [], 'ofdm_csi', 'ce', 'ce_mode', 'dftdd', ...
        'n_dd', 2, 'W_ici', 8, 'otfs_csi', 'ce', 'thr', 3, 'E', 10, 'tap_tol', 0.02, ...
        'mp_iter', 20, 'mp_damp', 0.7, 'mp_mode', 'edge', 'mp_tol', 1e-4);
    fn = fieldnames(d);
    for i = 1:numel(fn), if ~isfield(rx, fn{i}), rx.(fn{i}) = d.(fn{i}); end, end
    N = g.N;  M = g.M;
    mdm = modem(rx.modtype);
    A = mdm.alphabet * rx.data_amp;

    % ---- draws: channel, noise, data (in this order for CRN pairing) --------
    cr = chan_realize(ch, g, rs);
    Z  = sqrt(N0/2) * (randn(rs, N, M) + 1j*randn(rs, N, M));
    bits_all = randi(rs, [0 1], mdm.k, N*M);
    Xd = reshape(mdm.map(bits_all), N, M) * rx.data_amp;

    % ---- pilot budget --------------------------------------------------------
    if isempty(rx.guard_extra), guard_extra = ch.l_tau; else, guard_extra = rx.guard_extra; end
    G = ch.l_tau + 1 + guard_extra;
    n_nondata = G*M;
    switch rx.pilot
        case 'matched',  x_p = sqrt(n_nondata * 10^(rx.pilot_boost_dB/10)) * rx.data_amp;
        case 'explicit', x_p = rx.x_p;
        case 'none',     x_p = 0;
        otherwise, error('ber_frame:pilot', 'unknown pilot mode %s', rx.pilot);
    end
    E_pil = x_p^2;

    switch lower(wf)
        case 'cp-ofdm'
            if strcmp(rx.pilot, 'none') || ~strcmp(rx.ofdm_csi, 'ce')
                X = Xd;  masks.data = true(N, M);  masks.pilot = false(N, M);
                if ~strcmp(rx.pilot, 'none')    % keep the same overhead for fairness even with ideal CSI
                    [X, masks] = ofdm_pilot_lattice(Xd, g, ch, n_nondata, E_pil);
                end
            else
                [X, masks] = ofdm_pilot_lattice(Xd, g, ch, n_nondata, E_pil);
            end
            Y = io_cp_ofdm(X, g, cr, rx.use_ici) + Z;
            switch rx.ofdm_csi
                case 'ce'
                    H_hat = ofdm_ce(Y, X, masks, ch.l_tau, N0, mdm_scaled(mdm, rx.data_amp), rx.ce_mode, rx.n_dd);
                    Xs = ofdm_eq_1tap(Y, H_hat, N0);
                case 'ideal1tap'
                    % exact DIAGONAL of the (ICI) channel operator: per path the
                    % q = 0 kernel value I_p(0) times the no-ICI diagonal; ICI
                    % leakage is treated as noise by the one-tap MMSE.  (The
                    % no-ICI diagonal alone has the wrong phase/magnitude by
                    % I_p(0) - 1 = O(nu T) and floors the BER above dftdd CE.)
                    H_hat = zeros(N, M);
                    for q = 1:cr.P
                        crq = cr;  crq.alpha = cr.alpha(q);  crq.tau = cr.tau(q);  crq.nu = cr.nu(q);  crq.P = 1;
                        if rx.use_ici, I0 = ici_kernel(0, cr.nu(q), g.T); else, I0 = 1; end
                        H_hat = H_hat + I0 * io_cp_ofdm(ones(N, M), g, crq, false);
                    end
                    Xs = ofdm_eq_1tap(Y, H_hat, N0);
                case 'idealbanded'
                    H_hat = [];
                    Xs = ofdm_eq_banded(Y, g, cr, N0, rx.W_ici);
                otherwise, error('ber_frame:ofdm_csi', 'unknown ofdm_csi %s', rx.ofdm_csi);
            end
            out.H_hat = H_hat;
            taps = [];  mp_info = [];
        case 'cp-otfs'
            if strcmp(rx.pilot, 'none')
                X = Xd;  masks.data = true(N, M);  masks.known = false(N, M);  masks.known_val = zeros(N, M);
                masks.n_nondata = 0;
            else
                [X, masks] = otfs_pilot_frame(Xd, g, ch, x_p, guard_extra);
            end
            Y = io_cp_otfs(X, g, cr, rx.use_ici) + Z;
            switch rx.otfs_csi
                case 'ce'
                    assert(~strcmp(rx.pilot, 'none'), 'ber_frame:otfs_ce', 'OTFS CE needs a pilot');
                    taps = otfs_threshold_ce(Y, masks, x_p, N0, rx.thr);
                case 'ideal'
                    taps = dd_response_io(g, cr, rx.use_ici, rx.E, rx.tap_tol, 'origin');
                otherwise, error('ber_frame:otfs_csi', 'unknown otfs_csi %s', rx.otfs_csi);
            end
            [Xs, mp_info] = mp_detect_dd(Y, taps, A, N0, masks.known, masks.known_val, rx.mp_iter, rx.mp_damp, rx.mp_mode, rx.mp_tol);
            out.H_hat = [];
        otherwise
            error('ber_frame:wf', 'unknown waveform %s (BER receivers exist for cp-ofdm and cp-otfs)', wf);
    end

    % ---- hard decisions on the data cells -------------------------------------
    dm = masks.data;
    bits = mdm.demap(X(dm) / rx.data_amp);
    bits_hat = mdm.demap(Xs(dm) / rx.data_amp);
    Xhat = X;
    Xhat(dm) = reshape(mdm.map(bits_hat), [], 1) * rx.data_amp;
    out.bits = bits;  out.bits_hat = bits_hat;
    out.nbits = numel(bits);  out.nerr = nnz(bits ~= bits_hat);
    out.ber = out.nerr / out.nbits;
    out.ser = mean(Xhat(dm) ~= X(dm));
    out.X = X;  out.Xhat = Xhat;  out.Y = Y;  out.masks = masks;  out.cr = cr;
    out.taps = taps;  out.mp_info = mp_info;  out.x_p = x_p;  out.rx = rx;
end

function m = mdm_scaled(mdm, amp)
% modem whose map/demap operate at the scaled data amplitude
    m = mdm;
    m.alphabet = mdm.alphabet * amp;
    m.map   = @(bits) mdm.map(bits) * amp;
    m.demap = @(s) mdm.demap(s / amp);
end
