function idx = tau2idx(tau, Tsym, N)
    taus = getDelayArray(Tsym, N);
    taus = taus.';
    [~, idx] = min(abs(taus - tau));
end