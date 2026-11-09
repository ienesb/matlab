function Y = apply_dd_channel_eva(X, G)
% APPLY_DD_CHANNEL_EVA  Apply the effective DD response G to the DD symbol grid X
% as the circular 2D convolution
%   y[l,k] = sum_{l',k'} G[l',k'] x[[l-l']_N, [k-k']_M]        (noise added by caller).
% One circshift per nonzero tap -- exactly the forward model the MP detector
% assumes, so detection with these taps (build_effective_dd_eva) is model-matched.
    [N, M] = size(X);
    Y = zeros(N, M);
    [r, c, g] = find(G);
    for s = 1:numel(g)
        Y = Y + g(s) * circshift(X, [r(s)-1, c(s)-1]);
    end
end
