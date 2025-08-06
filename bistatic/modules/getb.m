function b = getb(tau, params)
    N = params.N;
    delta_f = params.delta_f;
    b = exp(-1j * 2*pi * (0:N-1)' * delta_f * tau);
end