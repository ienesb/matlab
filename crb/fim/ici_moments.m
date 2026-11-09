function mo = ici_moments(N, nu, g)
% ICI_MOMENTS  Data-independent ICI kernel sums of cp_ofdm_fim.tex eqs. 44-46
% and cp_ofdm_mcrb.tex (R):
%   S_k = sum_{n,n'} (n')^k |I(n'-n)|^2      k = 0,1,2   (real)
%   T_k = sum_{n,n'} (n')^k I*(n'-n) I'(n'-n) k = 0,1     (complex)
%   U   = sum_{n,n'} |I'(n'-n)|^2                        (real)
%   R   = sum_{n,n'} I*(n'-n) I''(n'-n)                  (complex; MCRB (nu,nu))
% Limits as I -> delta: S_k -> P_k, T_k -> 0, U -> 0, R -> 0.
    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    [I, Ip, Ipp] = ici_kernel(NP - NN, nu, g.T);   % rows n'
    np = (0:N-1).';
    rS = sum(abs(I).^2, 2);
    rT = sum(conj(I) .* Ip, 2);
    mo.S0 = sum(rS);  mo.S1 = sum(np .* rS);  mo.S2 = sum(np.^2 .* rS);
    mo.T0 = sum(rT);  mo.T1 = sum(np .* rT);
    mo.U  = sum(sum(abs(Ip).^2));
    mo.R  = sum(sum(conj(I) .* Ipp));
end
