function derivative = get_isi_coeffs_derivative(n, nu, T, tau)
% Doppler derivative of the partial kernel get_isi_coeffs
% (rcp_ofdm_fim.tex):
%   d/dnu I^<(n; nu, tau) = (1/T) int_0^tau j2pi t e^{j2pi (n Df + nu) t} dt
%                         = j2pi [ (tau/q) e^{q tau/T}
%                                  - (T/q^2) (e^{q tau/T} - 1) ]
% with q = j2pi (n + nu T).  At tau = T this reduces to the full-kernel
% derivative returned by get_ici_coeffs_derivative.
    q = 1j * 2*pi * (n + nu*T);
    E = exp(q * (tau/T));
    derivative = 1j*2*pi * ( (tau ./ q) .* E - (T ./ q.^2) .* (E - 1) );
    % Removable singularity at n + nu*T = 0: limit = j*pi*tau^2/T
    % (matches the j*pi*T limit of get_ici_coeffs_derivative at tau = T).
    derivative(q == 0) = 1j*pi*tau^2/T;
end
