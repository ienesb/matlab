function idx = tau2idx(tau, delta_f, N)
    taus = getDelayArray(delta_f, N);
    taus = taus.';
    [~, idx] = min(abs(taus - tau));
end