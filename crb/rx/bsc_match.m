function Xhat_b = bsc_match(X, Xhat, errmask, mdm, rs, mode)
% BSC_MATCH  Rate-matched, placement-randomized control arm: re-scatter the SAME
% number of symbol errors the real detector made, i.i.d. uniformly over the data
% cells.  Pilots/guard are kept exactly.  Isolates the error PLACEMENT from the
% error RATE (old/eva/mcrb_ce_ici.m bsc_match).
%   mode 'uniform' (default, legacy): each moved error goes to a uniformly random
%        DIFFERENT constellation point (for 4-QAM: E|e|^2 = 8/3, 1/3 antipodal);
%   mode 'nearest': flip ONE uniformly chosen Gray bit of the symbol (E|e|^2 = 2
%        for unit-power 4-QAM, the magnitude of a nearest-neighbour error, i.e.
%        the bsc_symbols('bit') model), so only the placement differs.
    if nargin < 5, rs = []; end
    if nargin < 6 || isempty(mode), mode = 'uniform'; end
    n_err = nnz(Xhat(errmask) ~= X(errmask));
    Xhat_b = X;
    if n_err == 0, return; end
    di = find(errmask);
    if isempty(rs), pe = di(randperm(numel(di), n_err)); else, pe = di(randperm(rs, numel(di), n_err)); end
    A = mdm.alphabet;
    switch lower(mode)
        case 'uniform'
            [~, lab] = min(abs(X(pe).' - A), [], 1);  lab = lab(:) - 1;
            if isempty(rs), r = randi(mdm.Q - 1, n_err, 1); else, r = randi(rs, mdm.Q - 1, n_err, 1); end
            Xhat_b(pe) = A(mod(lab + r, mdm.Q) + 1);
        case 'nearest'
            bits = mdm.demap(X(pe));                              % k x n_err (rows = bits)
            if isempty(rs), kb = randi(mdm.k, n_err, 1); else, kb = randi(rs, mdm.k, n_err, 1); end
            flip = false(size(bits));  flip(sub2ind(size(bits), kb(:), (1:n_err).')) = true;
            Xhat_b(pe) = mdm.map(xor(bits, flip));
        otherwise
            error('bsc_match:mode', 'unknown mode %s', mode);
    end
end
