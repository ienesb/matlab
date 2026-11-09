function e = relerr(A, B)
% RELERR  Relative Frobenius error ||A - B|| / max(||B||, eps).
    e = norm(A(:) - B(:)) / max(norm(B(:)), eps);
end
