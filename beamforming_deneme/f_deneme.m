Na = 64;
Ns = 1;
G = 90;

theta_min = -pi/4;
theta_max = pi/4;

n = (1:Na)';
phis = linspace(theta_min, theta_max, G);

A = exp(1j.*((n-1)*phis));

b = ones(G, 1);
% b(10:27) = 1;

% Given: A (matrix), b (vector)
f_unc = (A * A') \ (A * b); % Unconstrained solution

% Normalize to satisfy f^H A A^H f = 1
f_opt = f_unc / sqrt(f_unc' * (A * A') * f_unc);

f_opt'*A*A'*f_opt

norm(A'*f_opt-b)^2

figure;
plot(abs(A'*f_opt))
