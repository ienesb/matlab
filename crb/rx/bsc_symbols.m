function Xhat = bsc_symbols(X, p, mdm, rs, mode)
% BSC_SYMBOLS  Abstract detection-error model for the MCRB analysis.
%   Xhat = bsc_symbols(X, p, mdm, rs, mode)
%   mode 'bit' (default): each Gray label bit flips independently with
%        probability p (cp_ofdm_mcrb.tex BSC model; for Gray 4-QAM
%        E[Xhat|X] = (1-2p) X, E[E* Xhat] = -2p, E|E|^2 = 4p).
%   mode 'symbol': with probability p replace the symbol by a uniformly random
%        DIFFERENT constellation point (constellation-agnostic; used for 16-QAM).
% Entries not in the alphabet (pilots, zeros) are left untouched.
    if nargin < 4, rs = []; end
    if nargin < 5 || isempty(mode), mode = 'bit'; end
    sz = size(X);  x = X(:);
    A = mdm.alphabet;
    [dmin, idx] = min(abs(x.' - A), [], 1);
    isdata = dmin < 1e-9;
    Xhat = x;
    n = nnz(isdata);
    if n == 0, Xhat = reshape(Xhat, sz); return; end
    switch lower(mode)
        case 'bit'
            bits = mdm.demap(x(isdata));
            flip = urand(size(bits), rs) < p;
            Xhat(isdata) = mdm.map(xor(bits, flip));
        case 'symbol'
            err = urand([n 1], rs) < p;
            lab = idx(isdata).' - 1;
            newlab = mod(lab + irand(n, mdm.Q - 1, rs), mdm.Q);   % uniformly different
            lab(err) = newlab(err);
            Xhat(isdata) = A(lab + 1);
        otherwise
            error('bsc_symbols:mode', 'unknown mode %s', mode);
    end
    Xhat = reshape(Xhat, sz);
end

function u = urand(sz, rs)
    if isempty(rs), u = rand(sz); else, u = rand(rs, sz); end
end

function k = irand(n, K, rs)
    if isempty(rs), k = randi(K, n, 1); else, k = randi(rs, K, n, 1); end
end
