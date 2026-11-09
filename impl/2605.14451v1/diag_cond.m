% diag_cond.m - diagnose conditioning of H, J_bar
clear; clc;
N=128; L=40;
rng(20260514);
tau  = sort(127*rand(L,1));
fprintf('min gap = %.4f, max gap = %.4f, condH_a = %.2e\n', ...
    min(diff(tau)), max(diff(tau)), 0);

n_idx     = (0:N-1).';
a_mat     = exp(-1j*2*pi*(n_idx)*(tau.'/N));
adot_mat  = (-1j*2*pi*n_idx/N) .* a_mat;
beta = ones(L,1);
H = [adot_mat .* beta.', a_mat, 1j*a_mat];

fprintf('size(H) = %d x %d\n', size(H,1), size(H,2));
fprintf('rank(H) (numerical, complex)  = %d\n', rank(H));
fprintf('cond(H_a) = %.2e\n', cond(a_mat));
fprintf('cond(H) (real-augmented) = %.2e\n', cond([real(H);imag(H)]));

ReHH = real(H'*H);
fprintf('cond(Re{H^H H}) = %.2e\n', cond(ReHH));
fprintf('min eig (Re{H^H H}) = %.4e\n', min(eig(ReHH)));
