function [I, Ip, Ipp] = ici_kernel(q, nu, T)
% ICI_KERNEL  1/T-normalized (peak-1) OFDM ICI kernel and its REDUCED
% Doppler derivatives (latex/io/cp_ofdm_io.tex; latex/fim/cp_ofdm_fim.tex
% eqs. 2, 9; latex/mcrb/cp_ofdm_mcrb.tex line 84).
%
%   I(q)  = (1/T) int_0^T e^{j2pi (q df + nu) t} dt
%         = (e^{j2pi nu T} - 1) / (j2pi (q + nu T))
%   I'(q) := (1/j2pi) dI/dnu   = (1/T) int_0^T t   e^{...} dt
%   I''(q):= (1/j2pi) dI'/dnu  = (1/T) int_0^T t^2 e^{...} dt
%
% Inputs
%   q   integer lag n' - n (any array; NOT reduced mod N -- the kernel is
%       Toeplitz, not circulant).  For the N x N lag matrix use
%       [NP, NN] = ndgrid(0:N-1); q = NP - NN, so that Imat(n'+1, n+1) = I(n'-n)
%       and contractions read  Phi = Imat.' * A.
%   nu  Doppler [Hz] (scalar)
%   T   useful symbol duration 1/deltaf [s]
%
% Removable singularities at q + nu T = 0:  I = 1, I' = T/2, I'' = T^2/3.
% The full derivatives are dI/dnu = j2pi I', d^2I/dnu^2 = (j2pi)^2 I''.
% (The legacy get_ici_coeffs_derivative{,2} return the FULL derivatives.)
    D = 1j*2*pi*(q + nu*T);
    [m0, m1, m2] = exp_moments(D);
    I   = m0;
    Ip  = T   * m1;
    Ipp = T^2 * m2;
end
