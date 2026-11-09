function [A, B] = getAB_OFDM_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p)
% Single-target EXPECTED MCRB component matrices A and B for OFDM, ICI channel,
% i.i.d. unit-energy QPSK with BSC bit-error rate p (ofdm_mcrb.tex Sec. 1.1.2 /
% 2.1.2, "Expected ..." paragraphs).  Parameter order eta=[aR,aI,tau,nu].
% Depends on the true Doppler nu (through the kernel) but not on tau.
%
%   E{A} = -E{Jtilde} + E{A_Delta},     E{B} = E{Jtilde} + E{b} E{b}^T,
% with E{Jtilde} the expected ICI FIM and (eqs. 329-335, 210-213)
%   E{A_Delta}_{aR,tau} =  8pi p aI df M S1 /s2
%   E{A_Delta}_{aR,nu}  =  8pi p /s2 Im{alpha* W0}
%   E{A_Delta}_{aI,tau} = -8pi p aR df M S1 /s2
%   E{A_Delta}_{aI,nu}  =  8pi p /s2 Re{alpha* W0}
%   E{A_Delta}_{tau,tau}=  16pi^2 p |a|^2 M df^2 S2 /s2
%   E{A_Delta}_{nu,nu}  =  16pi^2 p |a|^2 V /s2
%   E{A_Delta}_{tau,nu} = -16pi^2 p |a|^2 df Re{W1} /s2
%   E{b} = [-4p M S0 aR/s2, -4p M S0 aI/s2, 0, 8pi p |a|^2 M Im{T0}/s2],
% using S_k, T_k, U, W0, W1, V (I' = (1/j2pi) dI/dnu convention).
%
% NOTE (documented gap): exact E{B} also needs Cov(b) (bias covariance), which
% ofdm_mcrb.tex does NOT provide.  This returns the mean-bias E{Jtilde}+E{b}E{b}^T.
%
% Inputs: grid (N,M), deltaf, Ts, sigma2, complex gain alpha, tau, nu, BSC p.
% Output: 4x4 A, B.

    aR = real(alpha); aI = imag(alpha); a2 = abs(alpha)^2;
    T  = 1/deltaf;

    % Kernel sums.  I' = (1/j2pi) dI/dnu , I'' = (1/j2pi) dI'/dnu.
    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    diff_n = NP - NN;
    I_mat   = get_ici_coeffs(diff_n, nu, T);
    Ip_mat  = get_ici_coeffs_derivative (diff_n, nu, T) / (1j*2*pi);
    Ipp_mat = get_ici_coeffs_derivative2(diff_n, nu, T) / (1j*2*pi)^2;
    npcol   = (0:N-1).';

    rowAbsI = sum(abs(I_mat).^2, 2);
    S0 = sum(rowAbsI);  S1 = sum(npcol.*rowAbsI);  S2 = sum(npcol.^2.*rowAbsI);
    rowTI = sum(conj(I_mat).*Ip_mat, 2);
    Tk0 = sum(rowTI);   Tk1 = sum(npcol.*rowTI);
    % R = sum I*(l) I''(l)  (from Theta = (mTs)^2 I + 2 mTs I' + I'').  The exact
    % expected Theta contraction uses Re{R}; ofdm_mcrb.tex eq.(317) writes
    % U = sum|I'|^2 instead, which differs from Re{R} at finite N (see note).
    R = sum(conj(I_mat(:)) .* Ipp_mat(:));

    M1 = M*(M-1)/2;
    M2 = M*(M-1)*(2*M-1)/6;
    W0 = Ts*M1*S0 + M*Tk0;
    W1 = Ts*M1*S1 + M*Tk1;
    V  = Ts^2*M2*S0 + Ts*M*(M-1)*real(Tk0) + M*real(R);

    % E{Jtilde} = expected ICI FIM.
    J = getFIM_OFDM_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, nu);

    % E{A_Delta}.
    Ad = zeros(4);
    Ad(1,3) =  8*pi*p*aI*deltaf*M*S1/sigma2;
    Ad(2,3) = -8*pi*p*aR*deltaf*M*S1/sigma2;
    Ad(1,4) =  8*pi*p/sigma2 * imag(conj(alpha)*W0);
    Ad(2,4) =  8*pi*p/sigma2 * real(conj(alpha)*W0);
    Ad(3,3) =  16*pi^2*p*a2*M*deltaf^2*S2/sigma2;
    Ad(4,4) =  16*pi^2*p*a2*V/sigma2;
    Ad(3,4) = -16*pi^2*p*a2*deltaf*real(W1)/sigma2;
    Ad(3,1)=Ad(1,3); Ad(4,1)=Ad(1,4); Ad(3,2)=Ad(2,3); Ad(4,2)=Ad(2,4); Ad(4,3)=Ad(3,4);

    A = -J + Ad;

    % E{b} and mean-bias E{B}.
    b = [ -4*p*M*S0*aR/sigma2;
          -4*p*M*S0*aI/sigma2;
           0;
           8*pi*p*a2*M*imag(Tk0)/sigma2 ];
    % Exact E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b).
    Cb = getMCRB_covb_OFDM(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, true);
    B = J + b*b.' + Cb;
end
