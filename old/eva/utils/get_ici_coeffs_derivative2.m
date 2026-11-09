function d2I = get_ici_coeffs_derivative2(n, nu, T)
% Second derivative d^2 I / d nu^2 of the normalized (unit-peak) ICI kernel
%   I(n) = (e^{j2pi nu T} - 1) / (j2pi (n + nu T)).
% Same kernel/convention as get_ici_coeffs and get_ici_coeffs_derivative
% (which returns the full first derivative dI/dnu).
%
% Closed form (with c = j2pi T, D = j2pi(n + nu T), E = e^{j2pi nu T}):
%   d2I = c^2 ( E D^2 - 2 E D + 2(E-1) ) / D^3,
% and the L'Hopital limit at the removable singularity n + nu*T = 0 is
%   d2I = c^2/3 = -(4 pi^2 T^2)/3.

    c = 1j*2*pi*T;
    D = 1j*2*pi*(n + nu*T);
    E = exp(1j*2*pi*nu*T);

    d2I = c.^2 .* (E.*D.^2 - 2*E.*D + 2*(E - 1)) ./ (D.^3);

    % Removable singularity n + nu*T = 0 (integer n at nu=0, etc.).
    d2I(D == 0) = c.^2 / 3;
end
