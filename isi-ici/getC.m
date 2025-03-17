function C = getC(nu, T, N, M)
    cM = exp(1j*2*pi*T*nu*(0:M-1)).';
    cICI = exp(1j*2*pi*T*nu*(0:N-1)/N).';
    
    c_nu = kron(cM, cICI);
    
    C = diag((c_nu));
end