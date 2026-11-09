function F = fields_rcp_otfs(x, g, tau, nu, use_ici, strict)
% FIELDS_RCP_OTFS  Per-target fields of the SAMPLED-receiver RCP-OTFS model,
% latex/io/rcp_otfs_io.tex eq. RCP_IO_general and latex/fim/rcp_otfs_fim.tex
% (sampled branch):
%
%   mu[l',k'] = (alpha/NM) e^{j2pi nu Ts l'} sum_{l,k} x[l,k] A(l'-l)
%               e^{-j2pi w(l') k/M} B(k-k')
%
%   A(d) = sum_n e^{j2pi n (d/N - df tau)}          delay Dirichlet kernel
%   B(d) = sum_m e^{j2pi m (d/M + T nu)}            Doppler kernel with T, not TS
%   w(l') = 1[l' < l_p], l_p = tau/Ts               boundary-wrap indicator
%
% The sampled branch carries NO OFDM-side ICI kernel: intra-symbol Doppler is
% the intra-frame phase e^{j2pi nu Ts l'} (latex/fim/README.md audit item 3:
% carrying both would double-count).  use_ici = true keeps that phase (exact
% model); use_ici = false drops it (eq. RCP_IO_noICI) and with it the Ts l' Phi
% term of the Doppler derivative; the wrap phase and B' are always retained.
%
% Outputs (N x M):  F.Phi, F.Psi (A' kernel, n'-weighted), F.Obar = Ts l' Phi + T Omega
%   d mu/d alpha_R = Phi, d mu/d tau = -j2pi df alpha Psi, d mu/d nu = j2pi alpha Obar.
%
% Requires a sample-spaced delay tau = l_p Ts (the branch is derived only for
% that case); off-grid tau errors out -- use snap=true / MT40-int or the
% continuous branch ('rcp-otfs-cont' = RCP-OFDM at ISFFT(x)).
%
% Ported from old/eva/utils/getFIM_RCP_OTFS_multi.m, whose use_ici=false path
% is this exact sampled model; its use_ici=true path (ICI kernel AND phase) is
% the double-counting hybrid and is not reproduced.
    if nargin < 6, strict = true; end   % strict=false only for finite-difference tests
    [N, M] = size(x);
    T  = g.T;  Ts = g.Ts;
    lp_real = tau/Ts;
    lp = round(lp_real);
    if strict
        assert(abs(lp_real - lp) < 1e-6, 'fields_rcp_otfs:offgrid', ...
            'Sampled RCP-OTFS needs tau = l_p*Ts (tau/Ts = %.6g). Use snap=true or ''rcp-otfs-cont''.', lp_real);
    end

    ns     = (0:N-1).';
    ms     = (0:M-1).';
    lprime = (0:N-1).';

    % Delay side: [l',k] = sum_l x[l,k] A^{(a)}(l'-l) via FFT along delay.
    Xf   = fft(x, [], 1);                                  % sum_l x e^{-j2pi n' l/N}
    C0   = exp(-1j*2*pi*g.deltaf*ns*tau) .* Xf;            % e^{-j2pi n' df tau} Xf
    C1   = ns .* C0;
    UK   = N * ifft(C0, [], 1);                            % sum_l A(l'-l) x[l,k]
    UK1  = N * ifft(C1, [], 1);                            % sum_l A'(l'-l) x[l,k]

    % Boundary-wrap phase on rows l' < l_p.
    wr = lprime < lp;
    wraprow = exp(-1j*2*pi*ms/M).';                        % 1 x M
    UK(wr,:)  = UK(wr,:)  .* wraprow;
    UK1(wr,:) = UK1(wr,:) .* wraprow;

    % Doppler kernels B(k-k'), B'(k-k') as M x M circulants, [B]_{k,k'}.
    [KK, KKp] = ndgrid(0:M-1, 0:M-1);
    Dd    = mod(KK - KKp, M);
    Eph   = exp(1j*2*pi*(ms*(0:M-1))/M);                   % [m, delta]
    gm    = exp(1j*2*pi*ms*T*nu);
    bvec  = (gm.'         * Eph).';                        % B(delta)
    bpvec = ((ms .* gm).' * Eph).';                        % B'(delta)
    Bm    = bvec (Dd + 1);
    Bpm   = bpvec(Dd + 1);

    if use_ici
        gph = exp(1j*2*pi*nu*Ts*lprime);                   % intra-frame Doppler phase
    else
        gph = ones(N, 1);
    end
    nrm = gph / (N*M);

    F.Phi = nrm .* (UK  * Bm);
    F.Psi = nrm .* (UK1 * Bm);
    Om    = nrm .* (UK  * Bpm);                            % T Omega / T
    if use_ici
        F.Obar = Ts * (lprime .* F.Phi) + T * Om;
    else
        F.Obar = T * Om;
    end
end
