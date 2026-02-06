function [tau_hat, nu_hat, X_hat] = gradient_step(tau_hat, nu_hat, alpha_hat, X_hat, Y)

N = 400;
M = 60;

deltaf = 120000;
Tsym = 1/deltaf*1.07; % 8.9167e-06;

lr_tau = 5e-10; % 5e-10; % 5e-22;
lr_nu = 1e12; % 1e12; % 1;
lr_X = 5e9; % 5e9; % 0.1;

reg_X = 1e-11;

ns = (0:N-1)';
ms = (0:M-1);

H_hat = 0;
for k = 1:1
    b_tau_k = getb(tau_hat, deltaf, N);
    c_nu_k = getc(nu_hat, Tsym, M);
    H_hat_k = alpha_hat(k)*(b_tau_k * c_nu_k.');
    H_hat = H_hat + H_hat_k;
end

Y_hat = H_hat .* X_hat;

% 2-D Hann window
wn = hann(N,'periodic'); wm = hann(M,'periodic');
W2 = wn * wm.';    
W2 = ones(N,M);

perr = W2 .* (Y - Y_hat);

dJ_dtau = -sum(4*pi*ns*deltaf.*imag(X_hat .* H_hat_k .* conj(perr)), "all");
dJ_dnu = -sum(4*pi*ms*Tsym.*imag(X_hat .* H_hat_k .* conj(perr)), "all");

% dJ_dX = -2 * real(H_hat .* conj(perr));
dJ_dX = -2 * conj(H_hat) .* perr;

tau_hat = tau_hat - lr_tau * dJ_dtau;
nu_hat = nu_hat - lr_nu * dJ_dnu;
X_hat = X_hat - lr_X * (dJ_dX + reg_X.*X_hat);

end