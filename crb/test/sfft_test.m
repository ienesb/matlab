M = 8; N = 16;
Xdd = randn(M,N) + 1j*randn(M,N);

Xtf = isfft(Xdd);
Xdd2 = sfft(Xtf);

fprintf("relative error = %.3e\n", norm(Xdd2 - Xdd,'fro')/norm(Xdd,'fro'));