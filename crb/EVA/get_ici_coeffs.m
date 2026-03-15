function coeff = get_ici_coeffs(n, nu, T)
    coeff = (exp(1j*2*pi*nu*T) - 1) ./ (1j*2*pi*(n + nu*T));
end