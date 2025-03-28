rho = 0.022;
D = rho * Drad + (1-rho) * Dcom;
D1 = D
[V,e] = eig(D);
b_opt = V(:, end);

b_opt.' * Drad * conj(b_opt)
norm(b_opt)

[V,e] = eig(Drad);
b_opt = V(:, end);

b_opt.' * Drad * conj(b_opt)
norm(b_opt)

