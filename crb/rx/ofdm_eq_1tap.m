function Xsoft = ofdm_eq_1tap(Y, H, N0)
% OFDM_EQ_1TAP  Per-bin one-tap MMSE equalizer  conj(H) Y / (|H|^2 + N0).
% (latex/main.tex eq. lmmse omits the conjugate -- do not copy it.)
    Xsoft = conj(H) .* Y ./ (abs(H).^2 + N0);
end
