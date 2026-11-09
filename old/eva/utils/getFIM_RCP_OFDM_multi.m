function J = getFIM_RCP_OFDM_multi(X, deltaf, Ts, sigma2, alpha, tau, nu, use_ici)
% Multi-target RCP-OFDM Fisher information matrix (rcp_ofdm_fim.tex).
% --------------------------------------------------------------------------
% RCP-OFDM = reduced-CP OFDM: the M OFDM symbols are sent back to back and
% a SINGLE cyclic prefix is prepended to the whole frame (vs. CP-OFDM,
% which prefixes every symbol; getFIM_OFDM_multi_*).  The frame CP makes
% the delay cyclic over the whole frame, so receive window m sees the tail
% of symbol m-1 (cyclically, window 0 sees symbol M-1) for t' < tau_p:
%
%   Y[n,m] = sum_p alpha_p e^{j2pi nu_p m T}
%              sum_{n'} e^{-j2pi n' Df tau_p}
%                ( X[n',m] I^>_p(n'-n) + X[n',m-1 mod M] I^<_p(n'-n) )
%          + Z[n,m]
%
% where the full ICI kernel splits exactly, I_p = I^>_p + I^<_p, into the
% [tau_p, T] (current-symbol) and [0, tau_p] (previous-symbol) partial
% integrals (get_ici_coeffs / get_isi_coeffs).  Delay therefore causes
% data-dependent ISI even at nu = 0 -- structural to the reduced CP -- and
% the ISI split ratio carries extra delay information.
%
% Sanity limits: M = 1 reproduces getFIM_OFDM_multi_* exactly (the split
% recombines on the same symbol); per realization,
% getFIM_RCP_OTFS_multi(x) should equal this FIM at X = ISFFT(x)
% (unitarity of the Wigner+SFFT receiver).
%
% Parameter ordering:  eta = [a1R,a1I,tau1,nu1, a2R,a2I,tau2,nu2, ...].
%
% Inputs:
%   X        : N x M transmitted TF symbol grid
%   deltaf   : subcarrier spacing [Hz]
%   Ts       : OFDM symbol period [s].  Accepted only for interface
%              compatibility; RCP-OFDM uses the useful duration T = 1/deltaf
%              (single frame CP, no per-symbol CP).
%   sigma2   : noise variance
%   alpha    : P x 1 complex path gains
%   tau      : P x 1 path delays   [s], 0 <= tau_p < T
%   nu       : P x 1 path Dopplers [Hz]
%   use_ici  : true  -> full Doppler ICI kernels and their nu-derivatives
%              false -> kernels frozen at nu = 0 (delay split retained,
%                       Doppler enters only via the per-window phase)
%
% Output: J is 4P x 4P.

    [N, M] = size(X);
    P      = numel(alpha);
    T      = 1/deltaf;          % useful duration; Ts intentionally unused

    ns = (0:N-1).';
    ms = (0:M-1);

    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    diff_n   = NP - NN;

    G = zeros(N*M, 4*P);

    for p = 1:P
        ap = alpha(p);  tp = tau(p);  nup = nu(p);
        if use_ici, nuk = nup; else, nuk = 0; end

        % Split kernels I = I^> + I^< evaluated at lag d = n - n'.
        I_lt = get_isi_coeffs(diff_n, nuk, T, tp);          % previous symbol
        I_gt = get_ici_coeffs(diff_n, nuk, T) - I_lt;       % current symbol

        % Boundary kernel for the tau-derivative:
        %   dI^</dtau = -dI^>/dtau = (1/T) e^{j2pi (d Df + nu) tau}.
        Ebnd = exp(1j*2*pi*(diff_n*deltaf + nuk)*tp) / T;

        % Data premultiplied by the delay phase ONLY (the per-window Doppler
        % phase belongs to the receive window m, not to the data symbol, so
        % it must be applied after the cyclic previous-symbol shift).
        Ap      = X .* exp(-1j*2*pi*deltaf*ns*tp);
        Ap_prev = circshift(Ap, 1, 2);          % column m <- column m-1 mod M

        ph = exp(1j*2*pi*T*ms*nup);             % per-window Doppler phase

        % Mean field (gradient w.r.t. alphaR is Phi, w.r.t. alphaI is j*Phi).
        Phi = ph .* (I_gt.' * Ap + I_lt.' * Ap_prev);

        % d mu / d tau  (without alpha): delay phase term + boundary term.
        Dtau = ph .* ( -1j*2*pi*deltaf * (I_gt.' * (ns .* Ap) + I_lt.' * (ns .* Ap_prev)) ...
                       + Ebnd.' * (Ap_prev - Ap) );

        % d mu / d nu  (without alpha): window phase term (+ kernel terms).
        Dnu = (1j*2*pi*T*ms) .* Phi;
        if use_ici
            dI_lt = get_isi_coeffs_derivative(diff_n, nuk, T, tp);
            dI_gt = get_ici_coeffs_derivative(diff_n, nuk, T) - dI_lt;
            Dnu   = Dnu + ph .* (dI_gt.' * Ap + dI_lt.' * Ap_prev);
        end

        col = 4*(p-1);
        G(:, col+1) =      Phi(:);
        G(:, col+2) = 1j * Phi(:);
        G(:, col+3) = ap * Dtau(:);
        G(:, col+4) = ap * Dnu(:);
    end

    J = (2/sigma2) * real(G' * G);
end
