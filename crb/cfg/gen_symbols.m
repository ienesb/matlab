function [S, bits, mdm] = gen_symbols(N, M, modtype, rs)
% GEN_SYMBOLS  i.i.d. uniform unit-power Gray-labelled symbol grid.
%   [S, bits, mdm] = gen_symbols(N, M, modtype, rs)
% S: N x M symbols; bits: k x (N*M) (column j = bits of S(j), column-major);
% rs: RandStream or [] (global stream).  See modem.m.
    mdm = modem(modtype);
    if nargin < 4, rs = []; end
    if isempty(rs)
        bits = randi([0 1], mdm.k, N*M);
    else
        bits = randi(rs, [0 1], mdm.k, N*M);
    end
    S = reshape(mdm.map(bits), N, M);
end
