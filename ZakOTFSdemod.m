function y = ZakOTFSdemod(x, M, N, padlen)
    x = reshape(x, [M + padlen, N]);
    x = x(1:M,:);
    x = x(:);
    y = zak(x, M);

end