function B = getB(tau, df, N, M)
    bN = exp(-1j*2*pi*df*tau*(0:N-1)).';
    bISI = exp(-1j*2*pi*df*tau*(0:M-1)/M).';
    
    b_tau = kron(bN, bISI);
    
    B = diag((b_tau));
end