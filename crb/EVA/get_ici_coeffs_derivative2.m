function d2I = get_ici_coeffs_derivative2(n, nu, T)
    % Second derivative of I with respect to nu
    % nu: frequency offset (scalar or vector)
    % T: symbol period
    % n: subcarrier index/offset
    
    term_inner = exp(1j * 2 * pi * nu * T) .* (2 * 1j * pi^2 * (n + nu * T).^2 - 2 * pi * (n + nu * T) - 1j) + 1j;
    num = T^2 * term_inner;
    den = pi * (n + nu * T).^3;
    
    d2I = num ./ den;
end