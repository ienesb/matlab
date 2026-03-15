function derivative = get_ici_coeffs_derivative(n, nu, T)
    temp1 = (1j*2*pi*T) * exp(1j*2*pi*nu*T);
    temp2 = 1j*2*pi*(n + nu*T);
    derivative = (temp1 * temp2 - temp1 + (1j*2*pi*T)) ./ temp2.^2;
end