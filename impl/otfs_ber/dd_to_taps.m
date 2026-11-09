function taps = dd_to_taps(G, tol)
% DD_TO_TAPS  Convert a (sparse) 2D effective DD response G (N x M) into the
% tap list consumed by the MP detector.
%
%   tol (optional): drop taps with |G| < tol*max|G| (relative pruning of the
%   Dirichlet inter-Doppler-interference tails).  Default 0 (keep all
%   nonzero entries).  Used for the ideal detector to bound complexity; the
%   estimated channel is already pruned by the detection threshold.
%
% Returns struct with column vectors:
%   dk   : Doppler shift k' (0..N-1) of each tap
%   dl   : delay shift l'   (0..M-1) of each tap
%   gain : complex tap coefficient G[k',l']

if nargin < 2, tol = 0; end
if tol > 0
    G(abs(G) < tol*max(abs(G(:)))) = 0;
end
[krow, lcol] = find(G ~= 0);
idx = sub2ind(size(G), krow, lcol);
taps.dk   = krow - 1;
taps.dl   = lcol - 1;
taps.gain = G(idx);
end
