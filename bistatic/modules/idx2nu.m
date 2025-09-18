function nu = idx2nu(idx, Tsym, M, lambda)
    nus = getDopplerArray(Tsym, M, lambda);
    nu = nus(idx);
end