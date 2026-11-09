function derivative = get_ici_coeffs_derivative(n, nu, T)
% Local override of crb/EVA/get_ici_coeffs_derivative that additionally
% handles the singular limit at  n + nu*T = 0  (i.e. the diagonal entry of
% the lag matrix at exactly nu = 0).  The L'Hopital limit of
%   I'(n; nu, T) := d/dnu [ (e^{j2pi nu T} - 1) / (j2pi (n + nu*T)) ]
% at (n, nu) = (0, 0) is j*pi*T.
%
% For all other points where temp2 != 0 the result matches the original
% closed-form derivative.

    temp1 = (1j*2*pi*T) * exp(1j*2*pi*nu*T);
    temp2 = 1j*2*pi*(n + nu*T);
    derivative = (temp1 .* temp2 - temp1 + (1j*2*pi*T)) ./ (temp2.^2);

    % L'Hopital limit when temp2 == 0 (n + nu*T = 0):
    derivative(temp2 == 0) = 1j*pi*T;
end
