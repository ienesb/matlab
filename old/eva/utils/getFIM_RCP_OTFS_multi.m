function J = getFIM_RCP_OTFS_multi(x, deltaf, Ts, sigma2, alpha, tau, nu, use_ici)
% Multi-target RCP-OTFS Fisher information matrix (delay-Doppler domain).
% --------------------------------------------------------------------------
% RCP-OTFS = reduced-CP OTFS: the M OFDM-like symbols are sent back to back
% and a SINGLE cyclic prefix is prepended to the whole frame (vs. CP-OTFS,
% which prefixes every symbol).  Per rcp_otfs_fim.tex the RCP-OTFS mean is
% the CP-OTFS delay-Doppler mean (cp_otfs_fim.tex) augmented by two
% deterministic unit-modulus phases and with the per-symbol period T_S
% replaced by the useful duration T (no per-symbol CP):
%
%   y[l',k'] = (1/NM) sum_p alpha_p e^{j2pi nu_p Tsamp l'}              (global phase)
%                * sum_{l,k} x[l,k] K_p(l',l)
%                  e^{-j2pi w_p(l') k/M}                                (wrap phase)
%                  B_p(k-k')  +  noise
%
% with the wrap indicator w_p(l') = 1[l' < l_p], l_p = round(tau_p/Tsamp),
% the OFDM-side delay kernel K_p (= Dirichlet A_p when there is no ICI, or
% the two-index ICI kernel from get_ici_coeffs when use_ici), and the
% Doppler Dirichlet kernel B_p(d) = sum_m e^{j2pi m (d/M + T nu_p)}.
%
% Parameter ordering:  eta = [a1R,a1I,tau1,nu1, a2R,a2I,tau2,nu2, ...].
%
% Inputs:
%   x        : N x M delay-Doppler symbol grid (delay rows, Doppler cols)
%   deltaf   : subcarrier spacing [Hz]
%   Ts       : OFDM symbol period [s].  Accepted only for interface
%              compatibility with the OFDM/CP-OTFS FIMs; RCP-OTFS uses the
%              useful duration T = 1/deltaf (single frame CP, no per-symbol CP).
%   sigma2   : noise variance
%   alpha    : P x 1 complex path gains
%   tau      : P x 1 path delays   [s]
%   nu       : P x 1 path Dopplers [Hz]
%   use_ici  : true  -> retain the OFDM-side ICI kernel I(n'-n)
%              false -> orthogonal subcarriers, I(n'-n) = delta[n'-n]
%
% Output: J is 4P x 4P.
%
% Implementation note (per realization vs. expectation): like the CP-OTFS
% derivation, the delay-Doppler Dirichlet kernels leak across bins, so a
% single QPSK realization yields a random FIM; its Monte-Carlo average over
% i.i.d. unit-energy data converges to the expected RCP-OTFS FIM of
% rcp_otfs_fim.tex.  run_fim_avg.m performs that averaging.

    [N, M] = size(x);
    P      = numel(alpha);
    T      = 1/deltaf;          % useful OFDM symbol duration  (= N*Tsamp)
    Tsamp  = 1/(N*deltaf);      % intra-symbol sample period   (= T/N)

    ns     = (0:N-1).';         % subcarrier index n'
    ms     = (0:M-1).';         % OFDM-symbol index m
    lprime = (0:N-1).';         % output delay index l'

    % Delay-domain DFT of the data once: Xf[n',k] = sum_l x[l,k] e^{-j2pi n' l/N}.
    Xf = fft(x, [], 1);

    % p-independent helpers.
    if use_ici
        [NP, NN] = ndgrid(0:N-1, 0:N-1);
        diff_n   = NP - NN;
    end
    [KK, KKp] = ndgrid(0:M-1, 0:M-1);
    Ddopp     = mod(KK - KKp, M);             % B circulant index map
    Ephase    = exp(1j*2*pi*(ms*(0:M-1))/M);  % [m, delta] for the Doppler kernels
    wraprow   = exp(-1j*2*pi*ms/M).';         % 1 x M boundary-wrap phase

    G = zeros(N*M, 4*P);

    for p = 1:P
        ap = alpha(p);  tp = tau(p);  np = nu(p);
        lp = round(tp / Tsamp);               % integer delay tap for wrap indicator

        % ---- delay side: K, K^{(1)} (and K^{I'} with ICI) via FFT -----------
        Ap_tf = exp(-1j*2*pi*deltaf*ns*tp) .* Xf;   % e^{-j2pi n' df tp} Xf[n',k]
        if use_ici
            I_mat  = get_ici_coeffs           (diff_n, np, T);
            Ip_mat = get_ici_coeffs_derivative(diff_n, np, T);
            C0  = I_mat.'  * Ap_tf;             % K       (subcarrier n)
            C1  = I_mat.'  * (ns .* Ap_tf);     % K^{(1)}  (n'-weighted)
            CIp = Ip_mat.' * Ap_tf;             % K^{I'}   (ICI nu-derivative)
        else
            C0  = Ap_tf;                        % K  = A   (Dirichlet)
            C1  = ns .* Ap_tf;                  % K1 = A'
        end
        UK  = N * ifft(C0, [], 1);              % [l',k] = sum_l K(l',l) x[l,k]
        UK1 = N * ifft(C1, [], 1);
        if use_ici, UKIp = N * ifft(CIp, [], 1); end

        % ---- boundary-wrap phase: rows l' < l_p get column k *= e^{-j2pi k/M}
        wr = (lprime < lp);                     % N x 1 logical (all false if lp==0)
        VK  = UK;   VK (wr,:) = UK (wr,:) .* wraprow;
        VK1 = UK1;  VK1(wr,:) = UK1(wr,:) .* wraprow;
        if use_ici, VKIp = UKIp;  VKIp(wr,:) = UKIp(wr,:) .* wraprow;  end

        % ---- Doppler Dirichlet kernels B(k-k'), B'(k-k') as M x M circulants
        gm    = exp(1j*2*pi*ms*T*np);           % e^{j2pi m T nu}, M x 1
        bvec  = (gm.'         * Ephase).';      % B(delta),  delta = 0..M-1
        bpvec = ((ms .* gm).' * Ephase).';      % B'(delta)
        BB    = bvec (Ddopp + 1);
        BBp   = bpvec(Ddopp + 1);

        % ---- assemble fields with global phase gph[l'] = e^{j2pi nu Tsamp l'}
        gph = exp(1j*2*pi*np*Tsamp*lprime);     % N x 1
        nrm = gph / (N*M);

        Phi = nrm .* (VK  * BB);
        Psi = nrm .* (VK1 * BB);
        if use_ici
            LamCP = nrm .* ( T*(VK*BBp) + (VKIp*BB) );
        else
            LamCP = nrm .* ( T*(VK*BBp) );
        end
        Obar = Tsamp * (lprime .* Phi) + LamCP; % augmented Doppler field

        col = 4*(p-1);
        G(:, col+1) =                     Phi(:);
        G(:, col+2) = 1j               *  Phi(:);
        G(:, col+3) = -1j*2*pi*deltaf*ap * Psi(:);
        G(:, col+4) =  1j*2*pi*ap        * Obar(:);
    end

    J = (2/sigma2) * real(G' * G);
end
