function J = gradient_error(tau_hat, nu_hat, alpha_hat, X_hat, Y)

N = 400;
M = 60;

deltaf = 120000;
Tsym = 1/deltaf*1.07; % 8.9167e-06;

H_hat = 0;
for k = 1:1 % 1:K
    b_tau_k = getb(tau_hat, deltaf, N);
    c_nu_k = getc(nu_hat, Tsym, M);
    H_hat_k = alpha_hat(k)*(b_tau_k * c_nu_k.');
    H_hat = H_hat + H_hat_k;
end

Y_hat = H_hat .* X_hat; % X_hat == X;

% 2-D Hann window
wn = hann(N,'periodic'); wm = hann(M,'periodic');
W2 = wn * wm.';    
W2 = ones(N,M);

perr = W2 .* (Y - Y_hat);

J = sum(abs(perr(:)).^2);

end