function coeff = get_ici_coeffs(n, nu, T)
    nuT   = nu * T;
    denom = 1j * 2*pi * (n + nuT);
    coeff = (exp(1j*2*pi*nuT) - 1) ./ denom;
    % L'Hopital limit = 1 when n + nu*T = 0 (0/0 case, integer n)
    coeff(denom == 0) = 1;
end