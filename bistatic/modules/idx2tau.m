function tau = idx2tau(idx, Tsym, N)
    taus = getDelayArray(Tsym, N);
    tau = taus(idx);
end