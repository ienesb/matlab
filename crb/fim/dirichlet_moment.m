function F = dirichlet_moment(a, z, K)
% DIRICHLET_MOMENT  F_a(z, K) = sum_{k=0}^{K-1} k^a e^{j2pi k z}  (direct sum;
% the closed form is 0/0 at z = 0 mod 1 where F_0 = K).  cp_ofdm_fim.tex eq. 80.
    k = 0:K-1;
    F = sum((k.^a) .* exp(1j*2*pi*k*z));
end
