function r = test_transforms()
% Unitary ISFFT/SFFT (Phase 0.3).
    r = [];
    N = 16; M = 8;
    x = randn(N, M) + 1j*randn(N, M);
    X = isfft_dd(x);
    r = tc(r, 'sfft(isfft(x)) == x',   relerr(sfft_dd(X), x), 1e-13);
    r = tc(r, 'isfft(sfft(X)) == X',   relerr(isfft_dd(sfft_dd(X)), X), 1e-13);
    r = tc(r, 'Parseval',              abs(norm(X(:))^2 - norm(x(:))^2)/norm(x(:))^2, 1e-13);
    % explicit definition, cp_otfs_io.tex eq. (1)
    Xd = zeros(N, M);
    for n = 0:N-1, for m = 0:M-1
        [L, K] = ndgrid(0:N-1, 0:M-1);
        Xd(n+1, m+1) = sum(sum(x .* exp(1j*2*pi*(m*K/M - n*L/N)))) / sqrt(N*M);
    end, end
    r = tc(r, 'isfft matches the eq. (1) double sum', relerr(X, Xd), 1e-12);
    oracle_paths('add', 'old');  c = onCleanup(@() oracle_paths('remove', 'old'));
    r = tc(r, 'isfft_dd == old isfft_u', relerr(X, isfft_u(x)), 1e-15);
    r = tc(r, 'sfft_dd == old sfft_u',   relerr(sfft_dd(X), sfft_u(X)), 1e-15);
end
