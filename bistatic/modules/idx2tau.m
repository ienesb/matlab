function tau = idx2tau(idx, delta_f, N)
    taus = getDelayArray(delta_f, N);
    tau = taus(idx);
end