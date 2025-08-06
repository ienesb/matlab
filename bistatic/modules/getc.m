function c = getc(nu, params)
    M = params.M;
    Tsym = params.Tsym;
    c = exp(1j * 2*pi * (0:M-1)' * Tsym * nu);
end