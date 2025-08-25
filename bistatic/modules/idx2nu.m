function nu = idx2nu(idx, delta_f, M)
    nus = getDopplerArray(delta_f, M, 1);
    nu = nus(idx);
end