function [Iisi, Iisi_p, Ebnd] = isi_kernel(q, nu, T, tau)
% ISI_KERNEL  Split (reduced-CP) kernel of latex/io/rcp_ofdm_io.tex and its
% derivatives (latex/fim/rcp_ofdm_fim.tex).
%
%   I^ISI(q) = (1/T) int_0^tau e^{j2pi (q df + nu) t} dt
%            = (e^{j2pi (q df + nu) tau} - 1) / (j2pi (q + nu T))
%   I^S(q)   = I(q) - I^ISI(q)          (compute as the difference; exact)
%   Iisi_p   := (1/j2pi) dI^ISI/dnu = (1/T) int_0^tau t e^{...} dt   (reduced)
%   Ebnd     := dI^ISI/dtau = (1/T) e^{j2pi (q df + nu) tau} = -dI^S/dtau
%
% Removable singularities at q + nu T = 0:  I^ISI = tau/T,  Iisi_p = tau^2/(2T).
% (The legacy get_isi_coeffs_derivative returns the FULL derivative j*pi*tau^2/T.)
%
% Inputs: q integer lag array, nu [Hz], T = 1/deltaf [s], tau [s] with 0 <= tau < T.
    Dp = 1j*2*pi*(q + nu*T) * (tau/T);       % = j2pi (q df + nu) tau
    [m0, m1] = exp_moments(Dp);
    Iisi   = (tau/T)   * m0;
    Iisi_p = (tau^2/T) * m1;
    Ebnd   = exp(Dp) / T;
end
