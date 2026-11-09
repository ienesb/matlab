rng(0);
N   = 1e6;
rho = 0.1;                          % per-bit BSC error prob

% Random Gray-coded QPSK
bits = randi([0 1], N, 2);
tx   = ((1 - 2*bits(:,1)) + 1j*(1 - 2*bits(:,2))) / sqrt(2);

rx = bsc_qpsk(tx, rho);

% Decompose: correct / adjacent (one bit off) / diagonal (both bits off)
correct  = mean(rx == tx);
diagonal = mean(rx == -tx);
adjacent = 1 - correct - diagonal;

fprintf('rho = %.3f (per-bit error)\n', rho);
fprintf('  P(correct)  = %.4f   theory (1-rho)^2   = %.4f\n', correct,  (1-rho)^2);
fprintf('  P(adjacent) = %.4f   theory 2*rho*(1-rho)= %.4f\n', adjacent, 2*rho*(1-rho));
fprintf('  P(diagonal) = %.4f   theory rho^2       = %.4f\n', diagonal, rho^2);

mean(conj(tx) .* rx)

tx2 = isfft(tx);
rx2 = isfft(rx);

mean(conj(tx2) .* rx2)
