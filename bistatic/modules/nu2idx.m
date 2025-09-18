function idx = nu2idx(nu, Tsym, M, lambda)
    nus = getDopplerArray(Tsym, M, lambda);
    nus = nus.';
    [~, idx] = min(abs(nus - nu));
end