function coeff = get_isi_coeffs(n, nu, T, tau)
% Partial (split) ICI/ISI kernel for reduced-CP waveforms (rcp_ofdm_fim.tex):
%   I^<(n; nu, tau) = (1/T) int_0^tau e^{j2pi (n Df + nu) t} dt
%                   = (e^{j2pi (n Df + nu) tau} - 1) / (j2pi (n + nu T))
% Same 1/T normalization as get_ici_coeffs, which is the full-interval
% kernel: I(n) = I^<(n) + I^>(n), with I^> the [tau, T] complement.
    nuT   = nu * T;
    denom = 1j * 2*pi * (n + nuT);
    coeff = (exp(denom * (tau/T)) - 1) ./ denom;
    % Removable singularity at n + nu*T = 0: limit = tau/T.
    coeff(denom == 0) = tau / T;
end
