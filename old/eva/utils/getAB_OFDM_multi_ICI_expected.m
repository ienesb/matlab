function [A, B] = getAB_OFDM_multi_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p)
% Multi-target EXPECTED MCRB component matrices A and B for OFDM, ICI channel,
% i.i.d. unit-energy QPSK with BSC bit-error rate p (ofdm_mcrb.tex Sec. 1.2.2 /
% 2.2.2).  Parameter stacking eta=[eta_1;...;eta_P], eta_p=[aR_p,aI_p,tau_p,nu_p].
%
%   E{A} = -E{Jtilde} + E{A_Delta},   E{B} = E{Jtilde} + E{b}E{b}^T.
% E{Jtilde} is the expected ICI multi FIM.  A_Delta is block-diagonal; off-diagonal
% E{A} blocks are -E{Jtilde}_{p,q}.  The diagonal blocks and bias use the
% cross-target QPSK+BSC inner products (mismatch summed over source target r):
%   E{<Phi_{E,r}, Xtilde_q>} = -2p * (cross moment), conjugation on r, with
%   S_a^{rq}=sum (n')^a I_r* I_q e^{-j2pi n' df (tau_q-tau_r)},  T_a^{rq} (I_q'),
%   R^{rq}=sum I_r* I_q'' ,  F_b^{rq}=sum_m m^b e^{j2pi m Ts (nu_q-nu_r)}.
% I' = (1/j2pi) dI/dnu , I'' = (1/j2pi) dI'/dnu.
%
% NOTE 1 (doc gap): exact E{B} also needs Cov(b); this returns the mean-bias form.
% NOTE 2 (doc): eqs.(372)-(381) give the diagonal A_Delta blocks "by analogy"
% (the cross-target -2p inner products are not written out); they are derived here
% as above, reducing to the single-target case at P=1.  The Theta (nu,nu) term
% uses R = sum I_r* I_q'' directly, not U = sum|I'|^2 (see getAB_OFDM_ICI_expected).
%
% Inputs: grid (N,M), deltaf, Ts, sigma2, alpha/tau/nu (P x 1), BSC p.
% Output: 4P x 4P A, B.

    P = numel(alpha);
    T = 1/deltaf;
    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    diff_n = NP - NN;
    npcol  = (0:N-1).';

    Imat = cell(P,1);  Ipmat = cell(P,1);  Ippmat = cell(P,1);
    for t = 1:P
        Imat{t}   = get_ici_coeffs(diff_n, nu(t), T);
        Ipmat{t}  = get_ici_coeffs_derivative (diff_n, nu(t), T) / (1j*2*pi);
        Ippmat{t} = get_ici_coeffs_derivative2(diff_n, nu(t), T) / (1j*2*pi)^2;
    end

    J = getFIM_OFDM_multi_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu);

    A = -J;                          % off-diagonal blocks final
    b = zeros(4*P, 1);
    c2 = 2/sigma2;  c4 = 4*pi/sigma2;  c8 = 8*pi^2/sigma2;

    for q = 1:P
        % Accumulate IP*_q = sum_r alpha_r^* E{<Phi_{E,r}, ._q>}.
        IPphi=0; IPpsi=0; IPxi=0; IPom=0; IPom1=0; IPth=0;
        for r = 1:P
            wph = exp(-1j*2*pi*deltaf*npcol*(tau(q)-tau(r)));     % N x 1
            rS = sum(conj(Imat{r}) .* Imat{q},   2);
            rT = sum(conj(Imat{r}) .* Ipmat{q},  2);
            rR = sum(conj(Imat{r}) .* Ippmat{q}, 2);
            S0 = sum(wph.*rS);  S1 = sum(npcol.*wph.*rS);  S2 = sum(npcol.^2.*wph.*rS);
            T0 = sum(wph.*rT);  T1 = sum(npcol.*wph.*rT);
            R0 = sum(wph.*rR);
            zM = Ts*(nu(q)-nu(r));
            F0 = dirichlet_moment(0, zM, M);
            F1 = dirichlet_moment(1, zM, M);
            F2 = dirichlet_moment(2, zM, M);

            w = conj(alpha(r)) * (-2*p);
            IPphi = IPphi + w *  S0*F0;
            IPpsi = IPpsi + w *  deltaf*S1*F0;
            IPxi  = IPxi  + w *  deltaf^2*S2*F0;
            IPom  = IPom  + w * (Ts*S0*F1 + T0*F0);
            IPom1 = IPom1 + w *  deltaf*(Ts*S1*F1 + T1*F0);
            IPth  = IPth  + w * (Ts^2*S0*F2 + 2*Ts*T0*F1 + R0*F0);
        end
        aq = alpha(q);

        Ad = zeros(4);
        Ad(1,3) =  c4*imag(IPpsi);   Ad(1,4) = -c4*imag(IPom);
        Ad(2,3) =  c4*real(IPpsi);   Ad(2,4) = -c4*real(IPom);
        Ad(3,3) = -c8*real(aq*IPxi);
        Ad(4,4) = -c8*real(aq*IPth);
        Ad(3,4) =  c8*real(aq*IPom1);
        Ad(3,1)=Ad(1,3); Ad(4,1)=Ad(1,4); Ad(3,2)=Ad(2,3); Ad(4,2)=Ad(2,4); Ad(4,3)=Ad(3,4);

        idx = 4*(q-1)+(1:4);
        A(idx, idx) = A(idx, idx) + Ad;

        b(idx) = [ c2*real(IPphi); -c2*imag(IPphi); c4*imag(aq*IPpsi); -c4*imag(aq*IPom) ];
    end

    % Exact E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b).
    Cb = getMCRB_covb_OFDM(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, true);
    B = J + b*b.' + Cb;
end

function F = dirichlet_moment(a, z, K)
    k = 0:K-1;
    F = sum( (k.^a) .* exp(1j*2*pi*k*z) );
end
