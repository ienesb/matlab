function G = build_effective_dd_eva(alpha, delay_taps, nu, N, M, deltaf, Emax)
% BUILD_EFFECTIVE_DD_EVA  Exact effective delay-Doppler channel response G
% (Raviteja et al., WCNC 2018, eq. 13), in the delay = rows (0..N-1) /
% Doppler = cols (0..M-1) convention used by the eva BER/MCRB scripts.
%
% Each path i has integer delay tap l_i, complex gain alpha_i, and fractional
% Doppler tap k_i = nu_i * M / deltaf (= beta_i + gamma_i).  A fractional Doppler
% spreads across the Doppler columns via the Dirichlet kernel; Emax truncates
% that inter-Doppler-interference spread to |q| <= Emax (Inf = exact).
% Used both to APPLY the channel (apply_dd_channel_eva) and to give the MP
% detector its tap list.
    G = zeros(N, M);
    kcol = (0:M-1);                         % output Doppler indices (columns)
    for i = 1:numel(alpha)
        li    = delay_taps(i);
        kreal = nu(i) * M / deltaf;
        beta  = round(kreal);
        gamma = kreal - beta;               % fractional Doppler in (-0.5, 0.5]
        hp    = alpha(i) * exp(-1j*2*pi*kreal*li/(N*M));   % delay-Doppler twist
        if abs(gamma) < 1e-9
            kidx = mod(beta, M);            % integer Doppler -> single column
            G(li+1, kidx+1) = G(li+1, kidx+1) + hp;
        else
            q = mod(beta - kcol, M);
            q = mod(q + M/2, M) - M/2;      % map to [-M/2, M/2)
            a = -q - gamma;
            row = hp * (exp(1j*2*pi*a) - 1) ./ (M*(exp(1j*(2*pi/M)*a) - 1));
            if isfinite(Emax), row(abs(q) > Emax) = 0; end
            G(li+1, :) = G(li+1, :) + row;  % spread across Doppler columns
        end
    end
end
