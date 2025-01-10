function y = ZakOTFSmod(x, M, N, padlen)
    y = izak(x);
    y = [reshape(y, [M, N]); zeros(padlen, N)];
    y = y(:);
end