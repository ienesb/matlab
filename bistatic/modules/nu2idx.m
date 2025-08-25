function idx = nu2idx(nu, delta_f, M)
    nus = getDopplerArray(delta_f, M, 1);
    nus = nus.';
    [~, idx] = min(abs(nus - nu));
end